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
#define SBTPTR sbt.s+sbt.o

char *sbuf_read(char *pathname);  /* defined in util.h */
char *dup_str(char *s, char *e);  /* dj: defined in dparser's util.h */

extern D_ParserTables parser_tables_gram;

unsigned int found_jac = 0, ini = 0, found_print = 0;


typedef struct symtab {
  char *ss;			/* symbol string: all vars*/
  char *de;             /* symbol string: all Des*/
  int deo[MXSYM];        /* offest of des */
  int vo[MXSYM];	/* offset of symbols */
  int lh[MXSYM];	/* lhs symbols? =9 if a state var*/
  int ini[MXSYM];        /* initial variable assignment =2 if there are two assignments */
  int ini0[MXSYM];        /* state initial variable assignment =2 if there are two assignments */
  int di[MXDER];	/* ith of state vars */
  int nv;			/* nbr of symbols */
  int ix;                       /* ith of curr symbol */
  int id;                       /* ith of curr symbol */
  int fn;			/* curr symbol a fn?*/
  int nd;			/* nbr of dydt */
  int pos;
  int pos_de;
  int ini_i; // #ini
  int statei; // # states
  int li; // # lhs
  int pi; // # param
} symtab;
symtab tb;

typedef struct sbuf {
  char s[MXBUF];	/* curr print buffer */
  int o;			/* offset of print buffer */
} sbuf;
sbuf sb;			/* buffer w/ current parsed & translated line */
        			/* to be stored in a temp file */
sbuf sbt; 


char *extra_buf, *model_prefix, *md5;
#ifndef __STANDALONE__
char *out2;
#endif


static FILE *fpIO, *fpIO2;


/* new symbol? if no, find it's ith */
int new_or_ith(const char *s) {
  int i, len, len_s=strlen(s);

  if (tb.fn) return 0;
  if (!strcmp("t", s)) return 0;
  if (!strcmp("time", s)) return 0;
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
  if (!strcmp("time",value)){
    sprintf(SBPTR, " t");
    sprintf(SBTPTR, "t");
    sb.o += 2;
    sbt.o += 1;
  } else {
    sprintf(SBPTR, " %s", value);
    sprintf(SBTPTR, "%s", value);
    sb.o += strlen(value)+1;
    sbt.o += strlen(value);

  }
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

  if (!strcmp("(", name) ||
      !strcmp(")", name) ||
      !strcmp(",", name)
      ) {
    sprintf(SBPTR, "%s",name);
    sb.o++;
    sprintf(SBTPTR,"%s",name);
    sbt.o++;

  }
  
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
    sprintf(SBTPTR,"=");
    sbt.o++;
  }
  
  if (!strcmp("|",name)){
    sprintf(SBPTR," ||");
    sb.o += 3;
    sprintf(SBTPTR,"||");
    sbt.o += 2;
  }

  if (!strcmp("&",name)){
    sprintf(SBPTR," &&");
    sb.o += 3;
    sprintf(SBTPTR,"&&");
    sbt.o += 2;
  }

  if (!strcmp("<>",name) ||
      !strcmp("~=",name) ||
      !strcmp("/=",name) 
      ){
    sprintf(SBPTR," !=");
    sb.o += 3;
    sprintf(SBTPTR,"!=");
    sbt.o += 2;
  }
  
  free(value);
  
  depth++;
  if (nch != 0) {
    
    
    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, " pow(");
      sb.o += 5;
    }
    for (i = 0; i < nch; i++) {      
      if (!strcmp("derivative", name) && i< 2) continue;
      if (!strcmp("der_rhs", name)    && i< 2) continue;
      if (!strcmp("derivative", name) && i==3) continue;
      if (!strcmp("der_rhs", name)    && i==3) continue;
      if (!strcmp("derivative", name) && i==4) continue;
      
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
      if (!strcmp("ini0", name)     && i == 1) continue;

      /* if (!strcmp("decimalint",name)){ */
      /* 	// Make implicit double */
      /* 	sprintf(SBPTR,".0"); */
      /* 	sb.o += 2; */
      /* } */

      tb.fn = (!strcmp("function", name) && i==0) ? 1 : 0;
      D_ParseNode *xpn = d_get_child(pn,i);
      wprint_parsetree(pt, xpn, depth, fn, client_data);
      
      if (!strcmp("print_command",name)){
	found_print = 1;
	char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if  (!strncmp(v,"print",5)){
	  fprintf(fpIO,"full_print;\n");
	  fprintf(fpIO2,"print;\n");
	} else {
	  fprintf(fpIO, "%s;\n", v);
	  fprintf(fpIO2,"%s;\n", v);
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
	    sprintf(sbt.s,"ode0_printf(");
	    sbt.o = 12;
 	  } else if (!strncmp(v,"jac0",4)) {
	    sprintf(sb.s,"JAC0_Rprintf(");
	    sb.o = 12;
	    sprintf(sbt.s,"jac0_printf(");
            sbt.o = 12;
          } else if (!strncmp(v,"ode",3)){
	    sprintf(sb.s,"ODE_Rprintf(");
            sb.o = 11;
	    sprintf(sbt.s,"ode_printf(");
            sbt.o = 11;
	  } else if (!strncmp(v,"jac",3)){
	    sprintf(sb.s,"JAC_Rprintf(");
            sb.o = 11;
	    sprintf(sbt.s,"jac_printf(");
            sbt.o = 11;
	  } else if (!strncmp(v,"lhs",3)){
	    sprintf(sb.s,"LHS_Rprintf(");
            sb.o = 11;
	    sprintf(sbt.s,"lhs_printf(");
            sbt.o = 11;
	  } else {
	    sprintf(sb.s,"Rprintf(");
            sb.o = 7;
	    sprintf(sbt.s,"printf(");
            sbt.o = 6;
	  }
        }
	if (i == 2){
	  sprintf(SBPTR,"%s",v);
	  sb.o = strlen(sb.s);
	  sprintf(SBTPTR,"%s",v);
          sbt.o = strlen(sbt.s);
	}
	if (i == 4){
	  fprintf(fpIO,  "%s;\n", sb.s);
	  fprintf(fpIO2, "%s;\n", sbt.s);
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
	  sprintf(SBTPTR,"df(%s)/dy(",v);
	} else {
	  // New statment
	  sb.o = 0;
	  sbt.o = 0;
          sprintf(sb.s,"__PDStateVar__[__CMT_NUM_%s__*(__NROWPD__)+",v);
	  sprintf(sbt.s,"df(%s)/dy(",v);
        }
	sb.o = strlen(sb.s);
	sbt.o = strlen(sbt.s);
        free(v);
	continue;
      }
      if ((!strcmp("jac",name)  || !strcmp("jac_rhs",name) ||
	   !strcmp("dfdy",name) || !strcmp("dfdy_rhs",name)) && i == 4){
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	sprintf(SBPTR, "__CMT_NUM_%s__]",v);
	sb.o = strlen(sb.s);
	sprintf(SBTPTR, "%s)",v);
	sbt.o = strlen(sbt.s);
	if (strcmp("jac",name) == 0 ||
	    strcmp("dfdy",name) == 0){
	  sprintf(SBPTR ," = ");
	  sb.o += 3;
	  sprintf(SBTPTR ,"=");
          sbt.o += 1;
        }
        free(v);
	continue;
      }
      
      //inits
      if (!strcmp("selection_statement", name) && i==1) {
        sprintf(sb.s, "if (");
        sb.o = strlen(sb.s);
	sprintf(sbt.s, "if (");
        sbt.o = strlen(sbt.s);
        continue;
      }
      if (!strcmp("selection_statement", name) && i==3) {
        sprintf(SBPTR, " {");
        sb.o += 2;
	sprintf(SBTPTR, "{");
        sbt.o += 1;
        fprintf(fpIO,  "%s\n", sb.s);
	fprintf(fpIO2, "%s\n", sbt.s);
        continue;
      }
      if (!strcmp("selection_statement__8", name) && i==0) {
        fprintf(fpIO,  "}\nelse {\n");
	fprintf(fpIO2, "}\nelse {\n");
        continue;
      }

      if (!strcmp("power_expression", name) && i==0) {
        sprintf(SBPTR, ",");
        sb.o++;
	sprintf(SBTPTR, "^");
        sbt.o++;
      }
      

      if (!strcmp("derivative", name) && i==2) {
        /* sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate[%d] +", tb.nd, tb.nd); */
        /* sb.o = strlen(sb.s); */
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if (new_de(v)){
	  sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate[%d] + ", tb.nd, tb.nd);
          sb.o = strlen(sb.s);
	  sprintf(sbt.s, "d/dt(%s)=", v);
	  sbt.o = strlen(sbt.s);
	  new_or_ith(v);
	  /* Rprintf("%s; tb.ini = %d; tb.ini0 = %d; tb.lh = %d\n",v,tb.ini[tb.ix],tb.ini0[tb.ix],tb.lh[tb.ix]); */
          if  ((tb.ini[tb.ix] == 1 && tb.ini0[tb.ix] == 0) || (tb.lh[tb.ix] == 1 & tb.ini[tb.ix] == 0)){
	    error("Cannot assign state variable %s; For initial condition assigment use '%s(0) = #'.\n",v,v);
	  }
          tb.lh[tb.ix] = 9;
          tb.di[tb.nd] = tb.ix;
          sprintf(tb.de+tb.pos_de, "%s,", v);
	  tb.pos_de += strlen(v)+1;
          tb.deo[++tb.nd] = tb.pos_de;
          /* free(buf); */
        } else {
	  new_or_ith(v);
          /* printf("de[%d]->%s[%d]\n",tb.id,v,tb.ix); */
          sprintf(sb.s, "__DDtStateVar__[%d] = ", tb.id);
	  sb.o = strlen(sb.s);
	  sprintf(sbt.s, "d/dt(%s)=", v);
          sbt.o = strlen(sbt.s);
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
	  sprintf(SBTPTR, "d/dt(%s)", v);
          sbt.o = strlen(sbt.s);
	}
        free(v);
	continue;
      }

      if ((!strcmp("assignment", name) || !strcmp("ini", name) || !strcmp("ini0", name)) && i==0) {
        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
	if (!strcmp("ini", name) || !strcmp("ini0", name)){
	  sprintf(sb.s, "(__0__)%s", v);
          sb.o = strlen(v)+7;
	  if (!strcmp("ini",name) & !new_de(v)){
	    error("Cannot assign state variable %s; For initial condition assigment use '%s(0) ='.\n",v,v);
	  }
        } else {
	  sprintf(sb.s, "%s", v);
          sb.o = strlen(v);
	  if (!new_de(v)){
	    error("Cannot assign state variable %s; For initial condition assigment use '%s(0) ='.\n",v,v);
	  }
        }
	sprintf(sbt.s, "%s", v);
        sbt.o = strlen(v);
        new_or_ith(v);
	if (!strcmp("assignment", name)){
	  tb.lh[tb.ix] = 1;
        } else if (!strcmp("ini", name) || !strcmp("ini0",name)){
	  if (tb.ini[tb.ix] == 0){
	    // If there is only one initialzation call, then assume
	    // this is a parameter with an initial value.
	    tb.ini[tb.ix] = 1;
	    if (!strcmp("ini0",name)){
	      tb.ini0[tb.ix] = 1;
	    }
	  } else {
	    // There is more than one call to this variable, it is a
	    // conditional variabile
	    tb.lh[tb.ix] = 1;
	    if (tb.ini0[tb.ix] == 1){
	      error("Cannot have conditional initial conditions for %s\n.",v);
	    }
	  }
	}
        free(v);
      }
    }

    if (!strcmp("assignment", name) || !strcmp("ini", name) || !strcmp("derivative", name) || !strcmp("jac",name) || !strcmp("dfdy",name) ||
	!strcmp("ini0",name)){
      fprintf(fpIO, "%s;\n", sb.s);
      fprintf(fpIO2, "%s;\n", sbt.s);
    }
    
    if (!strcmp("selection_statement", name)){
      fprintf(fpIO, "}\n");
      fprintf(fpIO2, "}\n");
    }
    
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
    retieve_var(i, buf);
    switch(scenario) {
    case 0:
      fprintf(outpt, i<tb.nv-1 ? "\t%s,\n" : "\t%s;\n", buf);
      break;
    case 1:
      fprintf(outpt, "\t%s = par_ptr[%d];\n", buf, j++);
      break;
    default: break;
    }
  }
  fprintf(outpt, "%s", post_str);  /* dj: security calls for const format */
}

void print_aux_info(FILE *outpt, char *model){
  int i, islhs,pi = 0,li = 0, o=0, statei = 0, ini_i = 0;
  char *s, *s2;
  char sLine[MXLEN+1];
  char buf[512], buf2[512];
  s = (char *) malloc(64*MXSYM);
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1){
      sprintf(s+o, "\tSET_STRING_ELT(lhs,%d,mkChar(\"%s\"));\n", li++, buf);
    } else if (strcmp(buf,"pi")){
      sprintf(s+o, "\tSET_STRING_ELT(params,%d,mkChar(\"%s\"));\n", pi++, buf);
    }
    o = strlen(s);
  }
  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    sprintf(s+o, "\tSET_STRING_ELT(state,%d,mkChar(\"%s\"));\n", statei++, buf);
    o = strlen(s);
  }
  fprintf(outpt,"extern SEXP %smodel_vars(){\n",model_prefix);
  fprintf(outpt,"\tSEXP lst    = PROTECT(allocVector(VECSXP, 7));\n");
  fprintf(outpt,"\tSEXP names  = PROTECT(allocVector(STRSXP, 7));\n");
  fprintf(outpt,"\tSEXP params = PROTECT(allocVector(STRSXP, %d));\n",pi);
  fprintf(outpt,"\tSEXP lhs    = PROTECT(allocVector(STRSXP, %d));\n",li);
  fprintf(outpt,"\tSEXP state  = PROTECT(allocVector(STRSXP, %d));\n",statei);
  fprintf(outpt,"\tSEXP tran   = PROTECT(allocVector(STRSXP, 7));\n");
  fprintf(outpt,"\tSEXP trann  = PROTECT(allocVector(STRSXP, 7));\n");
  fprintf(outpt,"\tSEXP mmd5   = PROTECT(allocVector(STRSXP, 2));\n");
  fprintf(outpt,"\tSEXP mmd5n  = PROTECT(allocVector(STRSXP, 2));\n");
  fprintf(outpt,"\tSEXP model  = PROTECT(allocVector(STRSXP, 3));\n");
  fprintf(outpt,"\tSEXP modeln = PROTECT(allocVector(STRSXP, 3));\n");
  fprintf(outpt,"%s",s);
  // Save for outputting in trans
  tb.pi = pi;
  tb.li = li;
  tb.statei = statei;
  fprintf(outpt,"\tSET_STRING_ELT(modeln,0,mkChar(\"model\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(model,0,mkChar(\"");
  for (i = 0; i < strlen(model); i++){
    if (model[i] == '"'){
      fprintf(outpt,"\\\"");
    } else if (model[i] == '\n'){
      fprintf(outpt,"\\n");
    } else if (model[i] == '\t'){
      fprintf(outpt,"\\t");
    } else if (model[i] >= 32  && model[i] <= 126){ // ASCII only
      fprintf(outpt,"%c",model[i]);
    }
  }
  fprintf(outpt,"\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(modeln,1,mkChar(\"normModel\"));\n");
  fpIO2 = fopen("out3.txt", "r");
  fprintf(outpt,"\tSET_STRING_ELT(model,1,mkChar(\"");
  err_msg((intptr_t) fpIO2, "Coudln't access out3.txt.\n", -1);
  while(fgets(sLine, MXLEN, fpIO2)) {  /* Prefered RxODE -- for igraph */
    for (i = 0; i < strlen(sLine); i++){
      if (sLine[i] == '"'){
        fprintf(outpt,"\\\"");
      } else if (sLine[i] == '\n'){
        fprintf(outpt,"\\n");
      } else if (sLine[i] == '\t'){
        fprintf(outpt,"\\t");
      } else if (sLine[i] >= 33  && sLine[i] <= 126){ // ASCII only
        fprintf(outpt,"%c",sLine[i]);
      }
    }
  }
  fclose(fpIO2);
  fprintf(outpt,"\"));\n");

#ifdef __STANDALONE__
  fpIO2 = fopen("out2.txt", "r");
#else
  fpIO2 = fopen(out2, "r");
#endif
  fprintf(outpt,"\tSET_STRING_ELT(modeln,2,mkChar(\"parseModel\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(model,2,mkChar(\"");
  err_msg((intptr_t) fpIO2, "Coudln't access out2.txt.\n", -1);
  while(fgets(sLine, MXLEN, fpIO2)) {  /* Prefered RxODE -- for igraph */
    for (i = 0; i < strlen(sLine); i++){
      if (sLine[i] == '"'){
        fprintf(outpt,"\\\"");
      } else if (sLine[i] == '\n'){
        fprintf(outpt,"\\n");
      } else if (sLine[i] == '\t'){
        fprintf(outpt,"\\t");
      } else if (sLine[i] >= 32  && sLine[i] <= 126){ // ASCII only
        fprintf(outpt,"%c",sLine[i]);
      }
    }
  }
  fclose(fpIO2);
  fprintf(outpt,"\"));\n");
#ifdef __STANDALONE__
  fpIO2 = fopen("out2.txt", "r");
#else
  fpIO2 = fopen(out2, "r");
#endif
  s[0] = '\0';
  o    = 0;
  while(fgets(sLine, MXLEN, fpIO2)) { 
    s2 = strstr(sLine,"(__0__)");
    if (s2){
      // See if this is a reclaimed initilization variable.
      for (i=0; i<tb.nv; i++) {
        if (tb.ini[i] == 1 && tb.lh[i] != 1){
          //(__0__)V2 =
          retieve_var(i, buf);
	  sprintf(buf2,"(__0__)%s =",buf);
          s2 = strstr(sLine,buf2);
          if (s2){
	    sprintf(s+o,"\tSET_STRING_ELT(inin,%d,mkChar(\"%s\"));\n",ini_i, buf);
	    o = strlen(s);
            sprintf(s+o,"\tREAL(ini)[%d] = %.*s;\n",ini_i++, strlen(sLine)-strlen(buf)-12,sLine + 10 + strlen(buf));
	    o = strlen(s);
            continue;
          }
        }
      }
      continue;
    }
  }
  fclose(fpIO2);
  // putin constants
  for (i=0; i<tb.nv; i++) {
    if (tb.ini[i] == 0 && tb.lh[i] != 1) {
      retieve_var(i, buf);
      // Put in constants
      if  (!strcmp("pi",buf)){
	sprintf(s+o,"\tSET_STRING_ELT(inin,%d,mkChar(\"pi\"));\n",ini_i);
	o = strlen(s);
	// Use well more digits than double supports
	sprintf(s+o,"\tREAL(ini)[%d] = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;\n",ini_i++);
	o = strlen(s);
      }
    }
  }
  tb.ini_i = ini_i;
  fprintf(outpt,"\tSEXP ini    = PROTECT(allocVector(REALSXP,%d));\n",ini_i);
  fprintf(outpt,"\tSEXP inin   = PROTECT(allocVector(STRSXP, %d));\n",ini_i);
  fprintf(outpt,"%s",s);
  // Vector Names
  fprintf(outpt,"\tSET_STRING_ELT(names,0,mkChar(\"params\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  0,params);\n");

  fprintf(outpt,"\tSET_STRING_ELT(names,1,mkChar(\"lhs\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  1,lhs);\n");
  
  fprintf(outpt,"\tSET_STRING_ELT(names,2,mkChar(\"state\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  2,state);\n");
  
  fprintf(outpt,"\tSET_STRING_ELT(names,3,mkChar(\"trans\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  3,tran);\n");
  
  fprintf(outpt,"\tSET_STRING_ELT(names,5,mkChar(\"model\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  5,model);\n");
  
  fprintf(outpt,"\tSET_STRING_ELT(names,4,mkChar(\"ini\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  4,ini);\n");

  fprintf(outpt,"\tSET_STRING_ELT(names,6,mkChar(\"md5\"));\n");
  fprintf(outpt,"\tSET_VECTOR_ELT(lst,  6,mmd5);\n");  
  
  // md5 values
  fprintf(outpt,"\tSET_STRING_ELT(mmd5n,0,mkChar(\"file_md5\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(mmd5,0,mkChar(\"%s\"));\n",md5);
  fprintf(outpt,"\tSET_STRING_ELT(mmd5n,1,mkChar(\"parsed_md5\"));\n");
  fprintf(outpt,"\tSET_STRING_ELT(mmd5,1,mkChar(__PARSED_MD5_STR__));\n");
  
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
  fprintf(outpt,"\tsetAttrib(model, R_NamesSymbol, modeln);\n");
  fprintf(outpt,"\tsetAttrib(ini, R_NamesSymbol, inin);\n");
  fprintf(outpt,"\tsetAttrib(lst, R_NamesSymbol, names);\n");

  fprintf(outpt,"\tUNPROTECT(13);\n");
  
  fprintf(outpt,"\treturn lst;\n");
  fprintf(outpt,"}\n");
  free(s);
  //fprintf(outpt,"SEXP __PARSED_MD5__()\n{\n\treturn %smodel_vars();\n}\n",model_prefix);
}

void codegen(FILE *outpt, int show_ode) {
  int i, j,print_ode=0, print_vars = 0, print_parm = 0, print_jac=0;
  char sLine[MXLEN+1];
  char buf[64];
  FILE *fpIO;

  char *hdft[]=
    {
      "#include <math.h>\n#ifdef __STANDALONE__\n#define Rprintf printf\n#define JAC_Rprintf printf\n#define JAC0_Rprintf if (jac_counter == 0) printf\n#define ODE_Rprintf printf\n#define ODE0_Rprintf if (dadt_counter == 0) printf\n#define LHS_Rprintf printf\n#define R_alloc calloc\n#else\n#include <R.h>\n#include <Rinternals.h>\n#include <Rmath.h>\n#define JAC_Rprintf Rprintf\n#define JAC0_Rprintf if (jac_counter == 0) Rprintf\n#define ODE_Rprintf Rprintf\n#define ODE0_Rprintf if (dadt_counter == 0) Rprintf\n#define LHS_Rprintf Rprintf\n#endif\n#define max(a,b) (((a)>(b))?(a):(b))\n#define min(a,b) (((a)<(b))?(a):(b))\nvoid update_par_ptr(double t);\n",
      "extern long dadt_counter;\nextern long jac_counter;\nextern double InfusionRate[99];\nextern double *par_ptr;\nextern double podo;\nextern double tlast;\n\n// prj-specific differential eqns\nvoid ",
      "dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n",
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
    fprintf(outpt,"\n");
    for (i = 0; i < strlen(extra_buf); i++){
      if (extra_buf[i] == '"'){
        fprintf(outpt,"\"");
      } else if (extra_buf[i] == '\n'){
        fprintf(outpt,"\n");
      } else if (extra_buf[i] == '\t'){
        fprintf(outpt,"\t");
      } else if (extra_buf[i] >= 32  && extra_buf[i] <= 126){ // ASCII only
        fprintf(outpt,"%c",extra_buf[i]);
      }
    }
    fprintf(outpt, "%s", hdft[1]);
    fprintf(outpt, "%s", model_prefix);
    fprintf(outpt, "%s", hdft[2]);
  } else if (show_ode == 2){
    fprintf(outpt, "// Jacobian derived vars\nvoid %scalc_jac(unsigned int neq, double t, double *__zzStateVar__, double *__PDStateVar__, unsigned int __NROWPD__) {\n\tdouble __DDtStateVar__[%d];\n",model_prefix,tb.nd+1);
  } else {
    fprintf(outpt, "// prj-specific derived vars\nvoid %scalc_lhs(double t, double *__zzStateVar__, double *lhs) {",model_prefix);
  }
  if (found_print){
    fprintf(outpt,"\n\tint __print_ode__ = 0, __print_vars__ = 0,__print_parm__ = 0,__print_jac__ = 0;\n");
  }
  if ((show_ode == 2 && found_jac == 1) || show_ode != 2){
    prnt_vars(0, outpt, 0, "double \n\t", "\n");     /* declare all used vars */
    fprintf(outpt,"\tupdate_par_ptr(t);\n");
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
      s = strstr(sLine,"(__0__)");
      if (s){
	// See if this is a reclaimed initilization variable.
	for (i=0; i<tb.nv; i++) {
	  if (tb.ini[i] == 1 && tb.lh[i] == 1){
	    //(__0__)V2 =
	    retieve_var(i, buf);
	    s = strstr(sLine,buf);
	    if (s){
	      fprintf(outpt,"\t%s\n",sLine + 7);
	      continue;
	    }
	  }
	}
	continue;
      }
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
    //fprintf(outpt,"\tfree(__ld_DDtStateVar__);\n");
    fprintf(outpt, "  jac_counter++;\n");
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
  memset(tb.lh,  0, MXSYM);
  memset(tb.ini, 0, MXSYM);
  memset(tb.ini0, 0, MXSYM);
  tb.nv=0;
  tb.nd=0;
  tb.fn=0;
  tb.ix=0;
  tb.id=0;
  tb.pos =0;
  tb.pos_de = 0;
  found_print = 0;
  found_jac = 0;
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
    fpIO  = fopen( "out2.txt", "w" );
#else
    fpIO = fopen( out2, "w" );
#endif
    fpIO2 = fopen( "out3.txt", "w" );
    err_msg((intptr_t) fpIO, "error opening out2.txt\n", -2);
    err_msg((intptr_t) fpIO2, "error opening out3.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    fclose(fpIO2);
    fpIO = fopen(c_file, "w");
    err_msg((intptr_t) fpIO, "error opening output c file\n", -2);
    codegen(fpIO, 1);
    codegen(fpIO, 2);
    codegen(fpIO, 0);
    print_aux_info(fpIO,buf);
    fclose(fpIO);
#ifdef __STANDALONE__
    remove("out2.txt");
    remove("out3.txt");
#endif
  } else {
    Rprintf("\nSyntax Error\n");
  }
#ifdef __STANDALONE__
  if (tb.ss) free(tb.ss);
  if (tb.de) free(tb.de);
  if (extra_buf) free(extra_buf);
  if (model_prefix) free(model_prefix);
#endif
}

#ifdef __STANDALONE__
int main(int argc, char *argv[]) {
  if (argc<3) {
    fprintf(stderr,"Usage: %s FILE_to_parse c_FILE [extra_c]\n",argv[0]);
    return -1;
  }
  model_prefix = (char *) malloc(1);
  model_prefix = '\0';
  printf("trans_internal(%s, %s)\n",argv[1],argv[2]);
  if (argc >= 3){ 
    extra_buf = sbuf_read(argv[3]); 
    if (!((intptr_t) extra_buf)){
      extra_buf = (char *) malloc(1);
      extra_buf[0] = '\0';
    }
  } else {
    if (!((intptr_t) extra_buf)){
      extra_buf = (char *) malloc(1);
      extra_buf[0] = '\0';
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
  const char *in, *out, *file, *pfile;
  char buf[512], buf2[512];
  char snum[512];
  char *s2;
  char sLine[MXLEN+1];
  int i, j, islhs, pi=0, li=0, ini_i = 0;
  double d;
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
      extra_buf = (char *) malloc(1);
      extra_buf[0]='\0';
    }
  } else {
    extra_buf = (char *) malloc(1); 
    extra_buf[0] = '\0';
  }
  
  if (model_prefix) free(model_prefix);
  if (isString(prefix) && length(prefix) == 1){
    model_prefix = CHAR(STRING_ELT(prefix,0));
  } else {
    model_prefix = (char *) malloc(1);
    model_prefix[0] = '\0';
  }

  if (md5) free(md5);
  if (isString(model_md5) && length(model_md5) == 1){
    md5 = CHAR(STRING_ELT(model_md5,0));
  } else {
    md5 = (char *) malloc(1);
    md5[0] = '\0';
  }
  
  if (out2) free(out2);
  if (isString(parse_model) && length(parse_model) == 1){
    out2 = CHAR(STRING_ELT(parse_model,0));
  } else {
    out2 = (char *) malloc(9); 
    sprintf(out2,"out2.txt"); 
  }
  trans_internal(in, out);
  SEXP lst   = PROTECT(allocVector(VECSXP, 6));
  SEXP names = PROTECT(allocVector(STRSXP, 6));
  
  SEXP tran  = PROTECT(allocVector(STRSXP, 7));
  SEXP trann = PROTECT(allocVector(STRSXP, 7));
  
  SEXP state = PROTECT(allocVector(STRSXP,tb.nd));
  
  SEXP params = PROTECT(allocVector(STRSXP, tb.pi));
  
  SEXP lhs    = PROTECT(allocVector(STRSXP, tb.li));
  
  SEXP inin   = PROTECT(allocVector(STRSXP, tb.ini_i));
  SEXP ini    = PROTECT(allocVector(REALSXP, tb.ini_i));
  
  SEXP model  = PROTECT(allocVector(STRSXP,3));
  SEXP modeln = PROTECT(allocVector(STRSXP,3));

  for (i=0; i<tb.nd; i++) {                     /* name state vars */
    retieve_var(tb.di[i], buf);
    SET_STRING_ELT(state,i,mkChar(buf));
  }
  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;      /* is a state var */
    retieve_var(i, buf);
    if (islhs == 1){
      SET_STRING_ELT(lhs,li++,mkChar(buf));
    } else if (strcmp(buf,"pi")){
      SET_STRING_ELT(params,pi++,mkChar(buf));
    } 
  }
  SET_STRING_ELT(names,0,mkChar("params"));
  SET_VECTOR_ELT(lst,  0,params);

  SET_STRING_ELT(names,1,mkChar("lhs"));
  SET_VECTOR_ELT(lst,  1,lhs);

  SET_STRING_ELT(names,2,mkChar("state"));
  SET_VECTOR_ELT(lst,  2,state);

  SET_STRING_ELT(names,3,mkChar("trans"));
  SET_VECTOR_ELT(lst,  3,tran);
  
  SET_STRING_ELT(names,4,mkChar("ini"));
  SET_VECTOR_ELT(lst,  4,ini);

  SET_STRING_ELT(names,5,mkChar("model"));
  SET_VECTOR_ELT(lst,  5,model);
  
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
  SET_STRING_ELT(tran,2,mkChar(buf)) ;

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

  fpIO2 = fopen(out2, "r");
  while(fgets(sLine, MXLEN, fpIO2)) { 
    s2 = strstr(sLine,"(__0__)");
    if (s2){
      // See if this is a reclaimed initilization variable.
      for (i=0; i<tb.nv; i++) {
        if (tb.ini[i] == 1 && tb.lh[i] != 1){
          //(__0__)V2 =
          retieve_var(i, buf);
	  sprintf(buf2,"(__0__)%s =",buf);
          s2 = strstr(sLine,buf2);
          if (s2){
	    /* Rprintf("%s[%d]->\n",buf,ini_i++); */
	    SET_STRING_ELT(inin,ini_i,mkChar(buf));
	    sprintf(snum,"%.*s",strlen(sLine)-strlen(buf)-12,sLine + 10 + strlen(buf));
	    sscanf(snum, "%lf", &d);
	    REAL(ini)[ini_i++] = d;
            continue;
          }
        }
      }
      continue;
    }
  }
  fclose(fpIO2);
  // putin constants
  for (i=0; i<tb.nv; i++) {
    if (tb.ini[i] == 0 && tb.lh[i] != 1) {
      retieve_var(i, buf);
      // Put in constants
      if  (!strcmp("pi",buf)){
	SET_STRING_ELT(inin,ini_i,mkChar("pi"));
	REAL(ini)[ini_i++] = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
      }
    }
  }
  file = sbuf_read(in);
  pfile = (char *) malloc(strlen(file)+1);
  j=0;
  for (i = 0; i < strlen(file); i++){
    if (file[i] == '"'  ||
	file[i] == '\n' ||
	file[i] == '\t' ||
	(file[i] >= 32 && file[i] <= 126)){
      sprintf(pfile+(j++),"%c",file[i]);
    }
  }
  SET_STRING_ELT(modeln,0,mkChar("model"));
  SET_STRING_ELT(model,0,mkChar(pfile));
  free(pfile);
  /* free(file); */
  SET_STRING_ELT(modeln,1,mkChar("normModel"));
  file = sbuf_read("out3.txt");
  if (file){
    pfile = (char *) malloc(strlen(file)+1);
    j=0;
    for (i = 0; i < strlen(file); i++){
      if (file[i] == '"'  ||
          file[i] == '\n' ||
          file[i] == '\t' ||
          (file[i] >= 33 && file[i] <= 126)){
        sprintf(pfile+(j++),"%c",file[i]);
      }
    }
    SET_STRING_ELT(model,1,mkChar(pfile));
    free(pfile);
  } else {
    SET_STRING_ELT(model,1,mkChar("Syntax Error"));
  }

  SET_STRING_ELT(modeln,2,mkChar("parseModel"));
  file = sbuf_read(out2);
  if (file){
    pfile = (char *) malloc(strlen(file)+1);
    j=0;
    for (i = 0; i < strlen(file); i++){
      if (file[i] == '"'  ||
          file[i] == '\n' ||
          file[i] == '\t' ||
          (file[i] >= 32 && file[i] <= 126)){
        sprintf(pfile+(j++),"%c",file[i]);
      }
    }
    SET_STRING_ELT(model,2,mkChar(pfile));
    free(pfile);
  } else {
    SET_STRING_ELT(model,2,mkChar("Syntax Error"));
  }
  
  
  setAttrib(ini,   R_NamesSymbol, inin);
  setAttrib(tran,  R_NamesSymbol, trann);
  setAttrib(lst,   R_NamesSymbol, names);
  setAttrib(model, R_NamesSymbol, modeln);
  UNPROTECT(11);
  remove("out3.txt");
  reset();
  return lst;
}


//FILE_to_parse c_FILE [aux_file_direcory extra_c]\n",argv[0]);
#endif
