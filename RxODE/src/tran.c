#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include "dparse_tree.h"
#ifdef __STANDALONE__
#define Rprintf printf
#define R_alloc calloc
#else
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#endif
#include "tran.g.d_parser.c"
#define max(a,b) (a)>(b) ? (a):(b)
#define MXSYM 5000
#define MXDER 500
#define MXLEN 1200
#define MXBUF 2400
#define SBPTR sb.s+sb.o

char *sbuf_read(char *pathname);  /* defined in util.h */
char *dup_str(char *s, char *e);  /* dj: defined in dparser's util.h */

extern D_ParserTables parser_tables_gram;

unsigned int found_jac = 0, ini = 0;


typedef struct symtab {
  char *ss;			/* symbol string: all vars*/
  char *de;             /* symbol string: all Des*/
  int deo[MXSYM];        /* offest of des */
  int vo[MXSYM];	/* offset of symbols */
  int lh[MXSYM];	/* lhs symbols? =9 if a state var*/
  int di[MXDER];	/* ith of state vars */
  int nv;			/* nbr of symbols */
  int ix;                       /* ith of curr symbol */
  int id;                       /* ith of curr symbol */
  int fn;			/* curr symbol a fn?*/
  int nd;			/* nbr of dydt */
  int pos;
  int pos_de;
} symtab;
symtab tb;

typedef struct sbuf {
  char s[MXBUF];	/* curr print buffer */
  int o;			/* offset of print buffer */
} sbuf;
sbuf sb;			/* buffer w/ current parsed & translated line */
        			/* to be stored in a temp file */

char *extra_buf, *model_prefix, *md5;
#ifndef __STANDALONE__
char *out2;
#endif


static FILE *fpIO;


/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=strlen(s);

  if (tb.fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("podo", s)) return 0;
  if (!strcmp("tlast", s)) return 0;
  if (!tb.nv) return 1;

  for (i=0; i<tb.nv; i++) {
    len = tb.vo[i+1] - tb.vo[i] - 1;  /* -1 for added ',' */
    if (!strncmp(tb.ss+tb.vo[i], s, max(len, len_s))) {	/* note we need take the max in order not to match a sub-string */
      tb.ix = i;
      return 0;
    }
  }
  return 1;
}

int new_de(const char *s){
  int i, len, len_s=strlen(s);
  for (i=0; i<tb.nd; i++) {
    len = tb.deo[i+1] - tb.deo[i] - 1;
    if (!strncmp(tb.de+tb.deo[i], s, max(len, len_s))) { /* note we need take the max in order not to match a sub-string */
      tb.id = i;
      return 0;
    }
  }
  return 1;
}

void wprint_node(int depth, char *name, char *value, void *client_data) {
  sprintf(SBPTR, " %s", value);
  sb.o += strlen(value)+1;
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  int nch = d_get_number_of_children(pn), i;
  char *value = (char*)dup_str(pn->start_loc.s, pn->end);
  char pexpr[80];


  if ((!strcmp("identifier", name) || !strcmp("identifier_no_output",name)) &&
      new_or_ith(value)) {
    /* printf("[%d]->%s\n",tb.nv,value); */
    sprintf(tb.ss+tb.pos, "%s,", value);
    tb.pos += strlen(value)+1;
    tb.vo[++tb.nv] = tb.pos;
    
  }

  if (!strcmp("(", name)) {sprintf(SBPTR, "("); sb.o++;}
  if (!strcmp(")", name)) {sprintf(SBPTR, ")"); sb.o++;}
  if (!strcmp(",", name)) {sprintf(SBPTR, ","); sb.o++;}
  
  if (
      !strcmp("identifier", name) ||
      !strcmp("constant", name) ||
      !strcmp("+", name) ||
      !strcmp("-", name) ||
      !strcmp("*", name) ||
      !strcmp("/", name) ||

      !strcmp("&&", name) ||
      !strcmp("||", name) ||
      !strcmp("!=", name) ||
      !strcmp("==", name) ||
      !strcmp("<=", name) ||
      !strcmp(">=", name) ||
      !strcmp("!", name) ||
      !strcmp("<", name) ||
      !strcmp(">", name) ||

      !strcmp("=", name)
     )
    fn(depth, name, value, client_data);

  // Operator synonyms  
  if (!strcmp("<-",name)){
    sprintf(SBPTR," =");
    sb.o += 2;
  }
  
  if (!strcmp("|",name)){
    sprintf(SBPTR," ||");
    sb.o += 3;
  }

  if (!strcmp("&",name)){
    sprintf(SBPTR," &&");
    sb.o += 3;
  }

  if (!strcmp("<>",name) ||
      !strcmp("~=",name) ||
      !strcmp("/=",name) 
      ){
    sprintf(SBPTR," !=");
    sb.o += 3;
  }
  
  free(value);
  
  depth++;
  if (nch != 0) {
    
    
    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, " pow(");
      sb.o += 5;
    }
    for (i = 0; i < nch; i++) {
      
      if (!strcmp("derivative", name) && i< 2)   continue;
      if (!strcmp("der_rhs", name)    && i< 2)   continue;
      if (!strcmp("derivative", name) && i==3)   continue;
      if (!strcmp("der_rhs", name)    && i==3) continue;
      if (!strcmp("derivative", name) && i==4)   continue;
      
      if (!strcmp("jac", name)     && i< 2)   continue;
      if (!strcmp("jac_rhs", name) && i< 2)   continue;
      if (!strcmp("jac", name)     && i == 3) continue;
      if (!strcmp("jac_rhs", name) && i == 3) continue;
      if (!strcmp("jac", name)     && i == 5) continue;
      if (!strcmp("jac_rhs", name) && i == 5) continue;
      if (!strcmp("jac", name)     && i == 6) continue;

      if (!strcmp("dfdy", name)     && i< 2)   continue;
      if (!strcmp("dfdy_rhs", name) && i< 2)   continue;
      if (!strcmp("dfdy", name)     && i == 3) continue;
      if (!strcmp("dfdy_rhs", name) && i == 3) continue;
      if (!strcmp("dfdy", name)     && i == 5) continue;
      if (!strcmp("dfdy_rhs", name) && i == 5) continue;
      if (!strcmp("dfdy", name)     && i == 6) continue;

      if (!strcmp("decimalint",name)){
	// Make implicit double
	sprintf(SBPTR,".0");
	sb.o += 2;
      }

      tb.fn = (!strcmp("function", name) && i==0) ? 1 : 0;
      D_ParseNode *xpn = d_get_child(pn,i);
      wprint_parsetree(pt, xpn, depth, fn, client_data);
      if (!strcmp("print_command",name)){
	char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if  (!strncmp(v,"print",5)){
	  fprintf(fpIO,"full_print;\n");
	} else {
	  fprintf(fpIO, "%s;\n", v);
	}
	/* sprintf(sb.s,"%s",v); */
        /* sb.o = str; */
        free(v);
      }
      if (!strcmp("printf_statement",name)){
	char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if (i == 0){
	  if (!strncmp(v,"ode0",4)){
	    sprintf(sb.s,"ODE0_Rprintf(");
            sb.o = 12;
 	  } else if (!strncmp(v,"jac0",4)) {
	    sprintf(sb.s,"JAC0_Rprintf(");
	    sb.o = 12;
          } else if (!strncmp(v,"ode",3)){
	    sprintf(sb.s,"ODE_Rprintf(");
            sb.o = 11;
	  } else if (!strncmp(v,"jac",3)){
	    sprintf(sb.s,"JAC_Rprintf(");
            sb.o = 11;
	  } else if (!strncmp(v,"lhs",3)){
	    sprintf(sb.s,"LHS_Rprintf(");
            sb.o = 11;
	  } else {
	    sprintf(sb.s,"Rprintf(");
            sb.o = 7;
	  }
        }
	if (i == 2){
	  sprintf(SBPTR,"%s",v);
	  sb.o = strlen(sb.s);
	}
	if (i == 4){
	  fprintf(fpIO, "%s;\n", sb.s);
  	}
	free(v);
	continue;
      } 

      if ( (!strcmp("jac",name) || !strcmp("jac_rhs",name) ||
	    !strcmp("dfdy",name) || !strcmp("dfdy_rhs",name)) && i == 2){
	found_jac = 1;
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if (!strcmp("jac_rhs",name) || !strcmp("dfdy_rhs",name)){
	  // Continuation statement
	  sprintf(SBPTR,"__PDStateVar__[__CMT_NUM_%s__*(__NROWPD__)+",v);
	} else {
	  // New statment
	  sb.o = 0;
          sprintf(sb.s,"__PDStateVar__[__CMT_NUM_%s__*(__NROWPD__)+",v);
        }
	sb.o = strlen(sb.s);
        free(v);
	continue;
      }
      if ((!strcmp("jac",name)  || !strcmp("jac_rhs",name) ||
	   !strcmp("dfdy",name) || !strcmp("dfdy_rhs",name)) && i == 4){
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	sprintf(SBPTR, "__CMT_NUM_%s__]",v);
	sb.o = strlen(sb.s);
	if (strcmp("jac",name) == 0){
	  sprintf(SBPTR ," = ");
	  sb.o += 3;
	}
        free(v);
	continue;
      }
      
      //inits
      if (!strcmp("selection_statement", name) && i==1) {
        sprintf(sb.s, "if (");
        sb.o = strlen(sb.s);
        continue;
      }
      if (!strcmp("selection_statement", name) && i==3) {
        sprintf(SBPTR, " {");
        sb.o += 2;
        fprintf(fpIO, "%s\n", sb.s);
        continue;
      }
      if (!strcmp("selection_statement__8", name) && i==0) {
        fprintf(fpIO, "}\nelse {\n");
        continue;
      }

      if (!strcmp("power_expression", name) && i==0) {
        sprintf(SBPTR, ",");
        sb.o++;
      }
      

      if (!strcmp("derivative", name) && i==2) {
        /* sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate[%d] +", tb.nd, tb.nd); */
        /* sb.o = strlen(sb.s); */
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if (new_de(v)){
	  sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate[%d] + ", tb.nd, tb.nd);
          sb.o = strlen(sb.s);
	  new_or_ith(v);
          tb.lh[tb.ix] = 9;
          tb.di[tb.nd] = tb.ix;
	  /* Rprintf("de[%d]->%s[%d]\n",tb.nd,v,tb.ix); */
          sprintf(tb.de+tb.pos_de, "%s,", v);
	  tb.pos_de += strlen(v)+1;
          tb.deo[++tb.nd] = tb.pos_de;
          /* free(buf); */
        } else {
	  new_or_ith(v);
          /* printf("de[%d]->%s[%d]\n",tb.id,v,tb.ix); */
          sprintf(sb.s, "__DDtStateVar__[%d] = ", tb.id);
	  sb.o = strlen(sb.s);
	}
        free(v);
        continue;
      }
      if (!strcmp("der_rhs", name)) {
      	char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        if (new_de(v)){
#ifdef __STANDALONE__
          fprintf(stderr,"Tried to use d/dt(%s) before it was defined",v);
	  free(v);
          if (tb.ss) free(tb.ss);
          if (tb.de) free(tb.de);
          exit(-1);
#else
	  error("Tried to use d/dt(%s) before it was defined",v);
	  free(v);
#endif
	} else {
	  sprintf(SBPTR, "__DDtStateVar__[%d]", tb.id);
	  sb.o = strlen(sb.s);
	}
        free(v);
	continue;
      }

      if (!strcmp("assignment", name) && i==0) {
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        sprintf(sb.s, "%s", v);
        sb.o = strlen(v);
        new_or_ith(v);
        tb.lh[tb.ix] = 1;
	free(v);
      }

    }

    if (!strcmp("assignment", name) || !strcmp("derivative", name) || !strcmp("jac",name) || !strcmp("dfdy",name))
      fprintf(fpIO, "%s;\n", sb.s);

    if (!strcmp("selection_statement", name))
      fprintf(fpIO, "}\n");

    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, ")");
      sb.o++;
    }
  }

}

void retieve_var(int i, char *buf) {
  int len;

  len = tb.vo[i+1] - tb.vo[i] - 1;
  strncpy(buf, tb.ss+tb.vo[i], len);
  buf[len] = 0;
}

void err_msg(int chk, const char *msg, int code)
{
  if(!chk) {
#ifdef __STANDALONE__
    if (tb.ss) free(tb.ss);
    if (tb.de) free(tb.de);
    fprintf(stderr, "%s", msg);
    exit(code);
#else
    error("%s",msg);
#endif
  }
}

/* when prnt_vars() is called, user defines the behavior in "case" */
void prnt_vars(int scenario, FILE *outpt, int lhs, const char *pre_str, const char *post_str) {
  int i, j;
  char buf[64];

  fprintf(outpt, "%s", pre_str);  /* dj: avoid security vulnerability */
  for (i=0, j=0; i<tb.nv; i++) {
    if (lhs && tb.lh[i]>0) continue;
    j++;
    retieve_var(i, buf);
    switch(scenario) {
    case 0: fprintf(outpt, i<tb.nv-1 ? "\t%s,\n" : "\t%s;\n", buf); break;
      case 1: fprintf(outpt, "\t%s = par_ptr[%d];\n", buf, j-1); break;
      default: break;
    }
  }
  fprintf(outpt, "%s", post_str);  /* dj: security calls for const format */
}

void print_aux_info(FILE *outpt){
  int i, islhs,pi = 0,li = 0, o=0, statei = 0;
  char *s;
  char buf[512];
  s = (char *) malloc(64*MXSYM);
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 0){
      sprintf(s+o, "\tSET_STRING_ELT(params,%d,mkChar(\"%s\"));\n", pi++, buf);
    } else {
      sprintf(s+o, "\tSET_STRING_ELT(lhs,%d,mkChar(\"%s\"));\n", li++, buf);
    }
    o = strlen(s);
  }
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    sprintf(s+o, "\tSET_STRING_ELT(state,%d,mkChar(\"%s\"));\n", statei++, buf);
    o = strlen(s);
  }
  fprintf(outpt,"extern SEXP %smodel_vars(){\n",model_prefix);
  fprintf(outpt,"\tSEXP lst = PROTECT(allocVector(VECSXP,5));\n");
  fprintf(outpt,"\tSEXP params = PROTECT(allocVector(STRSXP, %d));\n",pi);
  fprintf(outpt,"\tSEXP lhs    = PROTECT(allocVector(STRSXP, %d));\n",li);
  fprintf(outpt,"\tSEXP state  = PROTECT(allocVector(STRSXP, %d));\n",statei);
  fprintf(outpt,"\tSEXP names  = PROTECT(allocVector(STRSXP, 5));\n");
  fprintf(outpt,"\tSEXP tran   = PROTECT(allocVector(STRSXP, 7));\n");
  fprintf(outpt,"\tSEXP trann  = PROTECT(allocVector(STRSXP, 7));\n");
  fprintf(outpt,"\tSEXP mmd5   = PROTECT(allocVector(STRSXP, 2));\n");
  fprintf(outpt,"\tSEXP mmd5n  = PROTECT(allocVector(STRSXP, 2));\n");
  // Vector Names
  fprintf(outpt,"\tSET_STRING_ELT(names,0,mkChar(\"params\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,0,params);\n");
  fprintf(outpt,"\tSET_STRING_ELT(names,1,mkChar(\"lhs\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,1,lhs);\n");
  fprintf(outpt,"\tSET_STRING_ELT(names,2,mkChar(\"state\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,2,state);\n");
  fprintf(outpt,"\tSET_STRING_ELT(names,3,mkChar(\"trans\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,3,tran);\n");
  fprintf(outpt,"\tSET_STRING_ELT(names,4,mkChar(\"md5\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  4,mmd5);\n");

  // md5 values
  fprintf(outpt,"\tSET_STRING_ELT(mmd5n,0,mkChar(\"file_md5\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(mmd5,0,mkChar(\"%s\"));\n",md5);
  fprintf(outpt,"\tSET_STRING_ELT(mmd5n,1,mkChar(\"parsed_md5\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(mmd5,1,mkChar(__PARSED_MD5__));\n");

  fprintf(outpt,"%s",s);
  free(s);

  // now trans output
  fprintf(outpt,"\tSET_STRING_ELT(trann,0,mkChar(\"jac\"));\n");
  if (found_jac == 1){
    fprintf(outpt,"\tSET_STRING_ELT(tran,0,mkChar(\"fulluser\"));\n"); // Full User Matrix
  } else {
    fprintf(outpt,"\tSET_STRING_ELT(tran,0,mkChar(\"fullint\"));\n"); // Full Internal Matrix
  }
  fprintf(outpt,"\tSET_STRING_ELT(trann,1,mkChar(\"prefix\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 1,mkChar(\"%s\"));\n",model_prefix);

  fprintf(outpt,"\tSET_STRING_ELT(trann,2,mkChar(\"dydt\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 2,mkChar(\"%sdydt\"));\n",model_prefix);

  fprintf(outpt,"\tSET_STRING_ELT(trann,3,mkChar(\"calc_jac\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 3,mkChar(\"%scalc_jac\"));\n",model_prefix);

  fprintf(outpt,"\tSET_STRING_ELT(trann,4,mkChar(\"calc_lhs\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 4,mkChar(\"%scalc_lhs\"));\n",model_prefix);

  fprintf(outpt,"\tSET_STRING_ELT(trann,5,mkChar(\"model_vars\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 5,mkChar(\"%smodel_vars\"));\n",model_prefix);
  
  fprintf(outpt,"\tSET_STRING_ELT(trann,6,mkChar(\"ode_solver\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(tran, 6,mkChar(\"%sode_solver\"));\n",model_prefix);
  
  fprintf(outpt,"\tsetAttrib(tran, R_NamesSymbol, trann);\n");
  fprintf(outpt,"\tsetAttrib(mmd5, R_NamesSymbol, mmd5n);\n");
  fprintf(outpt,"\tsetAttrib(lst, R_NamesSymbol, names);\n");

  fprintf(outpt,"\tUNPROTECT(9);\n");
  
  fprintf(outpt,"\treturn lst;\n");
  fprintf(outpt,"}\n");
}

void codegen(FILE *outpt, int show_ode) {
  int i, j,print_ode=0, print_vars = 0, print_parm = 0, print_jac=0;
  char sLine[MXLEN+1];
  char buf[64];
  FILE *fpIO;

  char *hdft[]=
    {
      "#include <math.h>\n#ifdef __STANDALONE__\n#define Rprintf printf\n#define JAC_Rprintf printf\n#define JAC0_Rprintf if (jac_counter == 0) printf\n#define ODE_Rprintf printf\n#define ODE0_Rprintf if (dadt_counter == 0) printf\n#define LHS_Rprintf printf\n#define R_alloc calloc\n#else\n#include <R.h>\n#include <Rinternals.h>\n#include <Rmath.h>\n#define JAC_Rprintf Rprintf\n#define JAC0_Rprintf if (jac_counter == 0) Rprintf\n#define ODE_Rprintf Rprintf\n#define ODE0_Rprintf if (dadt_counter == 0) Rprintf\n#define LHS_Rprintf Rprintf\n#endif\n#define max(a,b) (((a)>(b))?(a):(b))\n#define min(a,b) (((a)<(b))?(a):(b))\n",
      "extern long dadt_counter;\nextern long jac_counter;\nextern double InfusionRate[99];\nextern double *par_ptr;\nextern double podo;\nextern double tlast;\n\n// prj-specific differential eqns\nvoid ",
      "dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n\tint __print_ode__ = 0, __print_vars__ = 0,__print_parm__ = 0,__print_jac__ = 0;\n",
      "    dadt_counter++;\n}\n\n"
    };
  if (show_ode == 1){
    fprintf(outpt, "%s", hdft[0]);
    if (found_jac == 1){
      for (i=0; i<tb.nd; i++) {                   /* name state vars */
        retieve_var(tb.di[i], buf);
        fprintf(outpt, "#define __CMT_NUM_%s__ %d\n", buf, i);
      }
    } 
    fprintf(outpt,"\n%s\n",extra_buf);
    fprintf(outpt, "%s", hdft[1]);
    fprintf(outpt, "%s", model_prefix);
    fprintf(outpt, "%s", hdft[2]);
  } else if (show_ode == 2){
    fprintf(outpt, "// Jacobian derived vars\nvoid %scalc_jac(unsigned int neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {\n\tint __print_ode__ = 0, __print_vars__ = 0,__print_parm__ = 0,__print_jac__ = 0;\n\tdouble __DDtStateVar__[%d];\n",model_prefix,tb.nd+1);
  } else {
    fprintf(outpt, "// prj-specific derived vars\nvoid %scalc_lhs(double t, double *__zzStateVar__, double *lhs) {\n\tint __print_ode__ = 0, __print_vars__ = 0,__print_parm__ = 0,__print_jac__ = 0;\n",model_prefix);
  }
  if ((show_ode == 2 && found_jac == 1) || show_ode != 2){
    prnt_vars(0, outpt, 0, "double", "\n");     /* declare all used vars */
    prnt_vars(1, outpt, 1, "", "\n");                   /* pass system pars */

    for (i=0; i<tb.nd; i++) {                   /* name state vars */
      retieve_var(tb.di[i], buf);
      fprintf(outpt, "\t%s = __zzStateVar__[%d];\n", buf, i);
    }
    fprintf(outpt,"\n");
#ifdef __STANDALONE__
    fpIO = fopen("out2.txt", "r");
#else
    fpIO = fopen(out2, "r");
#endif
    err_msg((intptr_t) fpIO, "Coudln't access out2.txt.\n", -1);
    while(fgets(sLine, MXLEN, fpIO)) {  /* parsed eqns */
      char *s;
      
      s = strstr(sLine,"ode_print;");
      if (show_ode == 1 && !s) s = strstr(sLine,"full_print;");
      if (show_ode != 1 && s) continue;
      else if (s) {
	fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
	fprintf(outpt,"\tRprintf(\"ODE Count: %%d\\tTime (t): %%f\\n\",dadt_counter,t);\n");
	fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
	fprintf(outpt,"\t__print_ode__ = 1;\n");
	fprintf(outpt,"\t__print_vars__ = 1;\n");
	fprintf(outpt,"\t__print_parm__ = 1;\n");
	print_ode  = 1;
	print_vars = 1;
	print_parm = 1;
	continue;
      }      
      s = strstr(sLine,"ODE_Rprintf");
      if ((show_ode != 1) && s) continue;
      
      s = strstr(sLine,"ODE0_Rprintf");
      if ((show_ode != 1) && s) continue;
      
      s = strstr(sLine, "__DDtStateVar__");
      if ((show_ode == 0) && s) continue;
      
      s = strstr(sLine,"JAC_Rprintf");
      if ((show_ode != 2) && s) continue;

      s = strstr(sLine,"jac_print;");
      if (show_ode == 2 && !s) s = strstr(sLine,"full_print;");
      if (show_ode != 2 && s) continue;
      else if (s) {
        fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
        fprintf(outpt,"\tRprintf(\"JAC Count: %%d\\tTime (t): %%f\\n\",jac_counter,t);\n");
        fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
	fprintf(outpt,"\t__print_ode__ = 1;\n");
	fprintf(outpt,"\t__print_jac__ = 1;\n");
        fprintf(outpt,"\t__print_vars__ = 1;\n");
        fprintf(outpt,"\t__print_parm__ = 1;\n");
        print_ode  = 1;
        print_vars = 1;
        print_parm = 1;
	print_jac = 1;
        continue;
      }

      s = strstr(sLine,"JAC0_Rprintf");
      if ((show_ode != 2) && s) continue;

      s = strstr(sLine,"LHS_Rprintf");
      if ((show_ode != 0) && s) continue;

      s = strstr(sLine,"lhs_print;");
      if (show_ode == 0 && !s) s = strstr(sLine,"full_print;");
      if (show_ode != 0 && s) continue;
      else if (s) {
        fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
        fprintf(outpt,"\tRprintf(\"LHS Time (t): %%f\\n\",t);\n");
        fprintf(outpt,"\tRprintf(\"================================================================================\\n\");\n");
        //fprintf(outpt,"\t__print_ode__ = 1;\n");
        fprintf(outpt,"\t__print_vars__ = 1;\n");
        fprintf(outpt,"\t__print_parm__ = 1;\n");
        //print_ode  = 1;
        print_vars = 1;
        print_parm = 1;
        continue;
      }
      
      s = strstr(sLine,"__PDStateVar__");
      if ((show_ode != 2) && s) continue;
      
      fprintf(outpt, "\t%s", sLine);
    }
    fclose(fpIO);
  }
  if (print_ode && show_ode != 0){
    fprintf(outpt,"\tif (__print_ode__ == 1){\n");
    for (i=0; i<tb.nd; i++) {                   /* name state vars */
      retieve_var(tb.di[i], buf);
      fprintf(outpt, "\t\tRprintf(\"d/dt(%s)[%d]:\\t%%f\\t%s:\\t%%f\\n\", __DDtStateVar__[%d],%s);\n", buf, i,buf,i,buf);
    }
    fprintf(outpt,"\t}\n");
  }
  if (print_jac && show_ode == 2){
    fprintf(outpt,"\tif (__print_jac__ == 1){\n");
    fprintf(outpt,"\tRprintf(\"Fixme\\n\");");
    fprintf(outpt,"\t}\n");
  }
  if (print_vars){
    fprintf(outpt,"\tif (__print_vars__ == 1){\n");
    fprintf(outpt,"\t\tRprintf(\".Left Handed Variables:.........................................................\\n\");\n");      
    for (i=0, j=0; i<tb.nv; i++) {
      if (tb.lh[i] != 1) continue;
      retieve_var(i, buf);
      fprintf(outpt, "\t\tRprintf(\"%s = %%f\\n\", %s);\n", buf, buf);
    }
    fprintf(outpt,"\t}\n");
  }
  if (print_parm){
    fprintf(outpt,"\tif (__print_parm__ == 1){\n");
    fprintf(outpt,"\t\tRprintf(\".User Supplied Variables:.......................................................\\n\");\n");
    for (i=0, j=0; i<tb.nv; i++) {
      if (tb.lh[i]>0) continue;
      j++;
      retieve_var(i, buf);
      fprintf(outpt, "\t\tRprintf(\"%s=%%f\\tpar_ptr[%d]=%%f\\n\",%s,par_ptr[%d]);\n", buf, j-1, buf,j-1);
    }
    fprintf(outpt,"\t}\n");
  }
  if (print_jac || print_vars || print_ode || print_parm){
    fprintf(outpt,"\tif (__print_jac__ || __print_vars__ || __print_ode__ || __print_parm__){\n");
    fprintf(outpt,"\t\tRprintf(\"================================================================================\\n\\n\\n\");\n\t}\n");
  }
  if (show_ode == 1){
    fprintf(outpt, "%s", hdft[3]);
  } else if (show_ode == 2){
    if (found_jac == 1){
      //fprintf(outpt,"\tfree(__ld_DDtStateVar__);\n");
      fprintf(outpt, "  jac_counter++;\n");
    }
    fprintf(outpt, "}\n");
  } else {
    fprintf(outpt, "\n");
    for (i=0, j=0; i<tb.nv; i++) {
      if (tb.lh[i] != 1) continue;
      retieve_var(i, buf);
      fprintf(outpt, "\tlhs[%d]=%s;\n", j, buf);
      j++;
    }
    fprintf(outpt, "}\n");
  }
}
void reset (){
  tb.vo[0]=0;
  tb.deo[0]=0;
  memset(tb.lh, 0, MXSYM);
  tb.nv=0;
  tb.nd=0;
  tb.fn=0;
  tb.ix=0;
  tb.id=0;
  tb.pos =0;
  tb.pos_de = 0;
}

void inits() {
  if (!ini){
    tb.ss = (char *) malloc(64*MXSYM);
    err_msg((intptr_t) tb.ss, "error allocating vars", 1);
    tb.de = (char *) malloc(64*MXSYM);
    err_msg((intptr_t) tb.de, "error allocating des", 1);
    ini = 1;
  }
  reset();
}


void trans_internal(char* parse_file, char* c_file){
  char *buf;
  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_gram, 1024);
  p->save_parse_tree = 1;
  buf = sbuf_read(parse_file);
  err_msg((intptr_t) buf, "error: empty buf for FILE_to_parse\n", -2);
  if ((pn=dparse(p, buf, strlen(buf))) && !p->syntax_errors) {
    inits();
#ifdef __STANDALONE__
    fpIO = fopen( "out2.txt", "w" );
#else
    fpIO = fopen( out2, "w" );
#endif
    err_msg((intptr_t) fpIO, "error opening out2.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    fpIO = fopen(c_file, "w");
    err_msg((intptr_t) fpIO, "error opening output c file\n", -2);
    codegen(fpIO, 1);
    codegen(fpIO, 2);
    codegen(fpIO, 0);
    print_aux_info(fpIO);
    fclose(fpIO);
#ifdef __STANDALONE__
    remove("out2.txt");
#endif
  } else {
    Rprintf("\nfailure\n");
  }
#ifdef __STANDALONE__
  if (tb.ss) free(tb.ss);
  if (tb.de) free(tb.de);
  if (extra_buf) free(extra_buf);
  if (model_prefix) free(model_prefix);
#else
  reset();
#endif
  
}

#ifdef __STANDALONE__
int main(int argc, char *argv[]) {
  if (argc<3) {
    fprintf(stderr,"Usage: %s FILE_to_parse c_FILE [extra_c]\n",argv[0]);
    return -1;
  }
  model_prefix = (char *) malloc(2);
  sprintf(model_prefix,"");
  printf("trans_internal(%s, %s)\n",argv[1],argv[2]);
  if (argc >= 3){ 
    extra_buf = sbuf_read(argv[3]); 
    if (!((intptr_t) extra_buf)){ 
      extra_buf = (char *) malloc(2); 
      sprintf(extra_buf,""); 
    }
  } else { 
    if (!((intptr_t) extra_buf)){ 
      extra_buf = (char *) malloc(2); 
      sprintf(extra_buf,""); 
    } 
  } 
     
  trans_internal(argv[1], argv[2]);
  return 0;
}

#else

void R_init_RxODE(DllInfo *info){
  inits();
}

void R_unload_RxODE(DllInfo *info){
  if (tb.ss) free(tb.ss);
  if (tb.de) free(tb.de);
  if (extra_buf) free(extra_buf);
  if (model_prefix) free(model_prefix);  
}

SEXP trans(SEXP parse_file, SEXP c_file, SEXP extra_c, SEXP prefix, SEXP model_md5,
	   SEXP parse_model){
  const char *in, *out;
  char buf[512];
  if (!isString(parse_file) || length(parse_file) != 1){
    error("parse_file is not a single string");
  }
  if (!isString(c_file) || length(c_file) != 1){
    error("c_file is not a single string");
  }
  in = CHAR(STRING_ELT(parse_file,0));
  out = CHAR(STRING_ELT(c_file,0));
  if (extra_buf) free(extra_buf);
  if (isString(extra_c) && length(extra_c) == 1){
    extra_buf = sbuf_read(CHAR(STRING_ELT(extra_c,0)));
    if (!((intptr_t) extra_buf)){ 
      extra_buf = (char *) malloc(2); 
      sprintf(extra_buf,""); 
    }
  } else {
    extra_buf = (char *) malloc(2); 
    sprintf(extra_buf,""); 
  }

  if (model_prefix) free(model_prefix);
  if (isString(prefix) && length(prefix) == 1){
    model_prefix = CHAR(STRING_ELT(prefix,0));
  } else {
    model_prefix = (char *) malloc(2);
    sprintf(model_prefix,"");
  }

  if (md5) free(md5);
  if (isString(model_md5) && length(model_md5) == 1){
    md5 = CHAR(STRING_ELT(model_md5,0));
  } else {
    md5 = (char *) malloc(2); 
    sprintf(md5,""); 
  }
  
  if (out2) free(out2);
  if (isString(parse_model) && length(parse_model) == 1){
    out2 = CHAR(STRING_ELT(parse_model,0));
  } else {
    out2 = (char *) malloc(9); 
    sprintf(out2,"out2.txt"); 
  }
  
  trans_internal(in, out);
  SEXP tran =PROTECT(allocVector(STRSXP, 8));
  SEXP trann=PROTECT(allocVector(STRSXP, 8));
  SET_STRING_ELT(trann,0,mkChar("jac"));
  if (found_jac == 1){
    SET_STRING_ELT(tran,0,mkChar("fulluser")); // Full User Matrix
  } else {
    SET_STRING_ELT(tran,0,mkChar("fullint")); // Full Internal Matrix
  }
  SET_STRING_ELT(trann,1,mkChar("prefix"));
  SET_STRING_ELT(tran,1,mkChar(model_prefix));

  sprintf(buf,"%sdydt",model_prefix);
  SET_STRING_ELT(trann,2,mkChar("dydt"));
  SET_STRING_ELT(tran,2,mkChar(buf));

  sprintf(buf,"%scalc_jac",model_prefix);
  SET_STRING_ELT(trann,3,mkChar("calc_jac"));
  SET_STRING_ELT(tran, 3,mkChar(buf));

  sprintf(buf,"%scalc_lhs",model_prefix);
  SET_STRING_ELT(trann,4,mkChar("calc_lhs"));
  SET_STRING_ELT(tran, 4,mkChar(buf));

  sprintf(buf,"%smodel_vars",model_prefix);
  SET_STRING_ELT(trann,5,mkChar("model_vars"));
  SET_STRING_ELT(tran, 5,mkChar(buf));

  sprintf(buf,"%sode_solver",model_prefix);
  SET_STRING_ELT(trann,6,mkChar("ode_solver"));
  SET_STRING_ELT(tran, 6,mkChar(buf));

  SET_STRING_ELT(trann,7,mkChar("file_md5"));
  SET_STRING_ELT(tran, 7,mkChar(md5));
  
  setAttrib(tran, R_NamesSymbol, trann);
  UNPROTECT(2);
  return tran;
}


//FILE_to_parse c_FILE [aux_file_direcory extra_c]\n",argv[0]);
#endif
