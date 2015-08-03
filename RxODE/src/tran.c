#include <stdio.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include "dparse_tree.h"
#define max(a,b) (a)>(b) ? (a):(b)
#define MXSYM 5000
#define MXDER 500
#define MXLEN 1200
#define MXBUF 2400
#define SBPTR sb.s+sb.o

char *sbuf_read(char *pathname);  /* defined in util.h */
char *dup_str(char *s, char *e);  /* dj: defined in dparser's util.h */
extern D_ParserTables parser_tables_gram;


typedef struct symtab {
  char *ss;			/* symbol string: all vars*/
  int vo[MXSYM];	/* offset of symbols */
  int lh[MXSYM];	/* lhs symbols? =9 if a state var*/
  int di[MXDER];	/* ith of state vars */
  int nv;			/* nbr of symbols */
  int ix;			/* ith of curr symbol */
  int fn;			/* curr symbol a fn?*/
  int nd;			/* nbr of dydt */
} symtab;
symtab tb;

typedef struct sbuf {
  char s[MXBUF];	/* curr print buffer */
  int o;			/* offset of print buffer */
} sbuf;
sbuf sb;			/* buffer w/ current parsed & translated line */
        			/* to be stored in a temp file */

static FILE *fpIO, *fp_inits;


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

void wprint_node(int depth, char *name, char *value, void *client_data) {
  sprintf(SBPTR, " %s", value);
  sb.o += strlen(value)+1;
}

void wprint_parsetree(D_ParserTables pt, D_ParseNode *pn, int depth, print_node_fn_t fn, void *client_data) {
  char *name = (char*)pt.symbols[pn->symbol].name;
  int nch = d_get_number_of_children(pn), i;
  char *value = (char*)dup_str(pn->start_loc.s, pn->end);
  char pexpr[80];


  if (!strcmp("identifier", name) && new_or_ith(value)) {
    static int pos=0;
    sprintf(tb.ss+pos, "%s,", value);
    pos += strlen(value)+1;
    tb.vo[++tb.nv] = pos;
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
  free(value);

  depth++;
  if (nch != 0) {

    if (!strcmp("power_expression", name)) {
      sprintf(SBPTR, " pow(");
      sb.o+=5;
    }

    for (i = 0; i < nch; i++) {
      if (!strcmp("derivative", name) && i< 2) continue;
      if (!strcmp("derivative", name) && i==3) continue;
      if (!strcmp("derivative", name) && i==4) continue;

      tb.fn = (!strcmp("function", name) && i==0) ? 1 : 0;
      D_ParseNode *xpn = d_get_child(pn,i);
      wprint_parsetree(pt, xpn, depth, fn, client_data);

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
        sprintf(sb.s, "__DDtStateVar__[%d] = InfusionRate[%d] +", tb.nd, tb.nd);
        sb.o = strlen(sb.s);

        char *v = (char*)dup_str(xpn->start_loc.s, xpn->end);
        new_or_ith(v);
        tb.lh[tb.ix] = 9;
        tb.di[tb.nd] = tb.ix;
        tb.nd++;
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

    if (!strcmp("assignment", name) || !strcmp("derivative", name))
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
    if (tb.ss) free(tb.ss);
    fprintf(stderr, "%s", msg);
    exit(code);
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

//output aux_files -- to be read by dvode() in R.
void prnt_aux_files(char *prefix) {
  int i, islhs;
  char buf[512];
  FILE *fp[3];

  sprintf(buf, "%sODE_PARS.txt",   prefix); fp[0] = fopen(buf, "w");
  sprintf(buf, "%sLHS_VARS.txt",   prefix); fp[1] = fopen(buf, "w");
  sprintf(buf, "%sSTATE_VARS.txt", prefix); fp[2] = fopen(buf, "w");
  i = (intptr_t) fp[0] * (intptr_t) fp[1] * (intptr_t) fp[2];  /* dj: sizeof(int)!=ptr */
  err_msg(i, "Coudln't open file to write.\n", -1);

  for (i=0; i<tb.nv; i++) {
    islhs = tb.lh[i];
    if (islhs>1) continue;	/* is a state var */
    retieve_var(i, buf);
    fprintf(fp[islhs], "%s ", buf);
  }

  for (i=0; i<tb.nd; i++) {			/* name state vars */
    retieve_var(tb.di[i], buf);
    fprintf(fp[2], "%s ", buf);
  }

  fclose(fp[0]);
  fclose(fp[1]);
  fclose(fp[2]);
}

void codegen(FILE *outpt, int show_ode) {
  int i, j;
  char sLine[MXLEN+1];
  char buf[64];
  FILE *fpIO;

  char *hdft[]=
    {
      "#include <math.h>\n#define max(a,b) (((a)>(b))?(a):(b))\n#define min(a,b) (((a)<(b))?(a):(b))\nextern long dadt_counter;\nextern double InfusionRate[99];\nextern double *par_ptr;\nextern double podo;\nextern double tlast;\n\n// prj-specific differential eqns\nvoid dydt(unsigned int neq, double t, double *__zzStateVar__, double *__DDtStateVar__)\n{\n",
      "    dadt_counter++;\n}\n\n"
    };

  if (show_ode)
    fprintf(outpt, "%s", hdft[0]);
  else
	fprintf(outpt, "// prj-specific derived vars\nvoid calc_lhs(double t, double *__zzStateVar__, double *lhs) {\n");
  prnt_vars(0, outpt, 0, "double\n", "\n");	/* declare all used vars */
  prnt_vars(1, outpt, 1, "", "\n");			/* pass system pars */

  for (i=0; i<tb.nd; i++) {			/* name state vars */
    retieve_var(tb.di[i], buf);
    fprintf(outpt, "\t%s = __zzStateVar__[%d];\n", buf, i);
  }
  fprintf(outpt,"\n");

  fpIO = fopen("out2.txt", "r");
  err_msg((intptr_t) fpIO, "Coudln't access out2.txt.\n", -1);
  while(fgets(sLine, MXLEN, fpIO)) {	/* parsed eqns */
	char *s;
	s = strstr(sLine, "InfusionRate");
	if (!show_ode && s) continue;
    fprintf(outpt, "\t%s", sLine);
  }

  if (show_ode) fprintf(outpt, "%s", hdft[1]);
  else {
  	fprintf(outpt, "\n");
    for (i=0, j=0; i<tb.nv; i++) {
      if (tb.lh[i] != 1) continue;
      retieve_var(i, buf);
      fprintf(outpt, "\tlhs[%d]=%s;\n", j, buf);
      j++;
    }
  	fprintf(outpt, "}\n");
  }

  fclose(fpIO);

}

void inits() {
  tb.ss = (char *) malloc(64*MXSYM);
  err_msg((intptr_t) tb.ss, "error allocating vars", 1);

  tb.vo[0]=0;
  memset(tb.lh, 0, MXSYM);
  tb.nv=0;
  tb.nd=0;
  tb.fn=0;
}

int main(int argc, char *argv[]) {
  char *buf;
  D_ParseNode *pn;
  /* any number greater than sizeof(D_ParseNode_User) will do;
     below 1024 is used */
  D_Parser *p = new_D_Parser(&parser_tables_gram, 1024);
  p->save_parse_tree = 1;

  if (argc<3) {
    fprintf(stderr,"Usage: %s FILE_to_parse c_FILE [aux file direcory]\n",argv[0]);
    return -1;
  } else {
    buf = sbuf_read(argv[1]);
    err_msg((intptr_t) buf, "error: empty buf\n", -2);
  }

  if ((pn=dparse(p, buf, strlen(buf))) && !p->syntax_errors) {
    inits();
    fpIO = fopen( "out2.txt", "w" );
    err_msg((intptr_t) fpIO, "error opening out2.txt\n", -2);
    wprint_parsetree(parser_tables_gram, pn, 0, wprint_node, NULL);
    fclose(fpIO);
    if (fp_inits) fclose(fp_inits);

	fpIO = fopen(argv[2], "w");
    codegen(fpIO, 1);
    codegen(fpIO, 0);
    fclose(fpIO);
    prnt_aux_files(argc<4 ? "" : argv[3]);
    remove("out2.txt");
    free(tb.ss);
  } else {
    printf("\nfailure\n");
  }
  return 0;
}

