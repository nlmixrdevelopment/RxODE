/*
Header file for using internal C-level dparser functions in RxODE (generated).
*/
#ifndef __RxODE_H__
#define __RxODE_H__
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Rdynload.h>
#include "gramgram.h"
#include "d.h"
#include "mkdparse.h"
#include "dparse.h"
#include "read_binary.h"
#if defined(__cplusplus)
extern "C" {
#endif
D_Parser *new_D_Parser(struct D_ParserTables *t, int sizeof_ParseNode_User){
  static D_Parser *(*fun)(struct D_ParserTables*, int)=NULL;
  if (fun == NULL) fun = (D_Parser* (*)(struct D_ParserTables*, int)) R_GetCCallable("RxODE","new_D_Parser");
  return fun(t, sizeof_ParseNode_User);
}

void   free_D_Parser (D_Parser *p ){
  static void   (*fun)(D_Parser* )=NULL;
  if (fun == NULL) fun = (void   (*)(D_Parser* )) R_GetCCallable("RxODE","free_D_Parser ");
  return fun(p );
}

D_ParseNode *dparse(D_Parser *p, char *buf, int buf_len){
  static D_ParseNode *(*fun)(D_Parser*, char*, int)=NULL;
  if (fun == NULL) fun = (D_ParseNode* (*)(D_Parser*, char*, int)) R_GetCCallable("RxODE","dparse");
  return fun(p, buf, buf_len);
}

void free_D_ParseNode(D_Parser *p, D_ParseNode *pn){
  static void (*fun)(D_Parser*, D_ParseNode*)=NULL;
  if (fun == NULL) fun = (void (*)(D_Parser*, D_ParseNode*)) R_GetCCallable("RxODE","free_D_ParseNode");
  return fun(p, pn);
}

void free_D_ParseTreeBelow(D_Parser *p, D_ParseNode *pn){
  static void (*fun)(D_Parser*, D_ParseNode*)=NULL;
  if (fun == NULL) fun = (void (*)(D_Parser*, D_ParseNode*)) R_GetCCallable("RxODE","free_D_ParseTreeBelow");
  return fun(p, pn);
}

int d_get_number_of_children(D_ParseNode *pn){
  static int (*fun)(D_ParseNode*)=NULL;
  if (fun == NULL) fun = (int (*)(D_ParseNode*)) R_GetCCallable("RxODE","d_get_number_of_children");
  return fun(pn);
}

D_ParseNode *d_get_child(D_ParseNode *pn, int child){
  static D_ParseNode *(*fun)(D_ParseNode*, int)=NULL;
  if (fun == NULL) fun = (D_ParseNode* (*)(D_ParseNode*, int)) R_GetCCallable("RxODE","d_get_child");
  return fun(pn, child);
}

D_ParseNode *d_find_in_tree(D_ParseNode *pn, int symbol){
  static D_ParseNode *(*fun)(D_ParseNode*, int)=NULL;
  if (fun == NULL) fun = (D_ParseNode* (*)(D_ParseNode*, int)) R_GetCCallable("RxODE","d_find_in_tree");
  return fun(pn, symbol);
}

char  * d_ws_before (D_Parser *p, D_ParseNode *pn ){
  static char  * (*fun)(D_Parser*, D_ParseNode* )=NULL;
  if (fun == NULL) fun = (char *  (*)(D_Parser*, D_ParseNode* )) R_GetCCallable("RxODE","d_ws_before ");
  return fun(p, pn );
}

char  * d_ws_after (D_Parser *p, D_ParseNode *pn ){
  static char  * (*fun)(D_Parser*, D_ParseNode* )=NULL;
  if (fun == NULL) fun = (char *  (*)(D_Parser*, D_ParseNode* )) R_GetCCallable("RxODE","d_ws_after ");
  return fun(p, pn );
}

void d_pass(D_Parser *p, D_ParseNode *pn, int pass_number){
  static void (*fun)(D_Parser*, D_ParseNode*, int)=NULL;
  if (fun == NULL) fun = (void (*)(D_Parser*, D_ParseNode*, int)) R_GetCCallable("RxODE","d_pass");
  return fun(p, pn, pass_number);
}

int resolve_amb_greedy(D_Parser *dp, int n, D_ParseNode **v){
  static int (*fun)(D_Parser*, int, D_ParseNode**)=NULL;
  if (fun == NULL) fun = (int (*)(D_Parser*, int, D_ParseNode**)) R_GetCCallable("RxODE","resolve_amb_greedy");
  return fun(dp, n, v);
}

char *d_dup_pathname_str(const char *str){
  static char *(*fun)(const char*)=NULL;
  if (fun == NULL) fun = (char* (*)(const char*)) R_GetCCallable("RxODE","d_dup_pathname_str");
  return fun(str);
}

void parse_whitespace(struct D_Parser *p, d_loc_t *loc, void **p_globals){
  static void (*fun)(struct D_Parser*, d_loc_t*, void**)=NULL;
  if (fun == NULL) fun = (void (*)(struct D_Parser*, d_loc_t*, void**)) R_GetCCallable("RxODE","parse_whitespace");
  return fun(p, loc, p_globals);
}

D_Scope *new_D_Scope(D_Scope *parent){
  static D_Scope *(*fun)(D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*)) R_GetCCallable("RxODE","new_D_Scope");
  return fun(parent);
}

D_Scope *enter_D_Scope(D_Scope *current, D_Scope *scope){
  static D_Scope *(*fun)(D_Scope*, D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*, D_Scope*)) R_GetCCallable("RxODE","enter_D_Scope");
  return fun(current, scope);
}

D_Scope *commit_D_Scope(D_Scope *scope){
  static D_Scope *(*fun)(D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*)) R_GetCCallable("RxODE","commit_D_Scope");
  return fun(scope);
}

D_Scope *equiv_D_Scope(D_Scope *scope){
  static D_Scope *(*fun)(D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*)) R_GetCCallable("RxODE","equiv_D_Scope");
  return fun(scope);
}

D_Scope *global_D_Scope(D_Scope *scope){
  static D_Scope *(*fun)(D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*)) R_GetCCallable("RxODE","global_D_Scope");
  return fun(scope);
}

D_Scope *scope_D_Scope(D_Scope *current, D_Scope *scope){
  static D_Scope *(*fun)(D_Scope*, D_Scope*)=NULL;
  if (fun == NULL) fun = (D_Scope* (*)(D_Scope*, D_Scope*)) R_GetCCallable("RxODE","scope_D_Scope");
  return fun(current, scope);
}

void free_D_Scope(D_Scope *st, int force){
  static void (*fun)(D_Scope*, int)=NULL;
  if (fun == NULL) fun = (void (*)(D_Scope*, int)) R_GetCCallable("RxODE","free_D_Scope");
  return fun(st, force);
}

D_Sym *new_D_Sym(D_Scope *st, char *name, char *end, int sizeof_D_Sym){
  static D_Sym *(*fun)(D_Scope*, char*, char*, int)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, char*, char*, int)) R_GetCCallable("RxODE","new_D_Sym");
  return fun(st, name, end, sizeof_D_Sym);
}

void free_D_Sym(D_Sym *sym){
  static void (*fun)(D_Sym*)=NULL;
  if (fun == NULL) fun = (void (*)(D_Sym*)) R_GetCCallable("RxODE","free_D_Sym");
  return fun(sym);
}

D_Sym *find_D_Sym(D_Scope *st, char *name, char *end){
  static D_Sym *(*fun)(D_Scope*, char*, char*)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, char*, char*)) R_GetCCallable("RxODE","find_D_Sym");
  return fun(st, name, end);
}

D_Sym *find_global_D_Sym(D_Scope *st, char *name, char *end){
  static D_Sym *(*fun)(D_Scope*, char*, char*)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, char*, char*)) R_GetCCallable("RxODE","find_global_D_Sym");
  return fun(st, name, end);
}

D_Sym *update_D_Sym(D_Sym *sym, D_Scope **st, int sizeof_D_Sym){
  static D_Sym *(*fun)(D_Sym*, D_Scope**, int)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Sym*, D_Scope**, int)) R_GetCCallable("RxODE","update_D_Sym");
  return fun(sym, st, sizeof_D_Sym);
}

D_Sym *update_additional_D_Sym(D_Scope *st, D_Sym *sym, int sizeof_D_Sym){
  static D_Sym *(*fun)(D_Scope*, D_Sym*, int)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, D_Sym*, int)) R_GetCCallable("RxODE","update_additional_D_Sym");
  return fun(st, sym, sizeof_D_Sym);
}

D_Sym *current_D_Sym(D_Scope *st, D_Sym *sym){
  static D_Sym *(*fun)(D_Scope*, D_Sym*)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, D_Sym*)) R_GetCCallable("RxODE","current_D_Sym");
  return fun(st, sym);
}

D_Sym *find_D_Sym_in_Scope(D_Scope *st, D_Scope *cur, char *name, char *end){
  static D_Sym *(*fun)(D_Scope*, D_Scope*, char*, char*)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope*, D_Scope*, char*, char*)) R_GetCCallable("RxODE","find_D_Sym_in_Scope");
  return fun(st, cur, name, end);
}

D_Sym *next_D_Sym_in_Scope(D_Scope **st, D_Sym **sym){
  static D_Sym *(*fun)(D_Scope**, D_Sym**)=NULL;
  if (fun == NULL) fun = (D_Sym* (*)(D_Scope**, D_Sym**)) R_GetCCallable("RxODE","next_D_Sym_in_Scope");
  return fun(st, sym);
}

void print_scope(D_Scope *st){
  static void (*fun)(D_Scope*)=NULL;
  if (fun == NULL) fun = (void (*)(D_Scope*)) R_GetCCallable("RxODE","print_scope");
  return fun(st);
}

Grammar *new_D_Grammar(char *pathname){
  static Grammar *(*fun)(char*)=NULL;
  if (fun == NULL) fun = (Grammar* (*)(char*)) R_GetCCallable("RxODE","new_D_Grammar");
  return fun(pathname);
}

void free_D_Grammar(Grammar *g){
  static void (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*)) R_GetCCallable("RxODE","free_D_Grammar");
  return fun(g);
}

int build_grammar(Grammar *g){
  static int (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*)) R_GetCCallable("RxODE","build_grammar");
  return fun(g);
}

int parse_grammar(Grammar *g, char *pathname, char *str){
  static int (*fun)(Grammar*, char*, char*)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*, char*, char*)) R_GetCCallable("RxODE","parse_grammar");
  return fun(g, pathname, str);
}

void print_grammar(Grammar *g){
  static void (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*)) R_GetCCallable("RxODE","print_grammar");
  return fun(g);
}

void print_rdebug_grammar(Grammar *g, char *pathname){
  static void (*fun)(Grammar*, char*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, char*)) R_GetCCallable("RxODE","print_rdebug_grammar");
  return fun(g, pathname);
}

void print_states(Grammar *g){
  static void (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*)) R_GetCCallable("RxODE","print_states");
  return fun(g);
}

void print_rule(Rule *r){
  static void (*fun)(Rule*)=NULL;
  if (fun == NULL) fun = (void (*)(Rule*)) R_GetCCallable("RxODE","print_rule");
  return fun(r);
}

void print_term(Term *t){
  static void (*fun)(Term*)=NULL;
  if (fun == NULL) fun = (void (*)(Term*)) R_GetCCallable("RxODE","print_term");
  return fun(t);
}

Production *lookup_production(Grammar *g, char *name, int len){
  static Production *(*fun)(Grammar*, char*, int)=NULL;
  if (fun == NULL) fun = (Production* (*)(Grammar*, char*, int)) R_GetCCallable("RxODE","lookup_production");
  return fun(g, name, len);
}

Rule *new_rule(Grammar *g, Production *p){
  static Rule *(*fun)(Grammar*, Production*)=NULL;
  if (fun == NULL) fun = (Rule* (*)(Grammar*, Production*)) R_GetCCallable("RxODE","new_rule");
  return fun(g, p);
}

Elem *new_elem_nterm(Production *p, Rule *r){
  static Elem *(*fun)(Production*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(Production*, Rule*)) R_GetCCallable("RxODE","new_elem_nterm");
  return fun(p, r);
}

void new_declaration(Grammar *g, Elem *e, uint kind){
  static void (*fun)(Grammar*, Elem*, uint)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, Elem*, uint)) R_GetCCallable("RxODE","new_declaration");
  return fun(g, e, kind);
}

Production *new_production(Grammar *g, char *name){
  static Production *(*fun)(Grammar*, char*)=NULL;
  if (fun == NULL) fun = (Production* (*)(Grammar*, char*)) R_GetCCallable("RxODE","new_production");
  return fun(g, name);
}

Elem *new_string(Grammar *g, char *s, char *e, Rule *r){
  static Elem *(*fun)(Grammar*, char*, char*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(Grammar*, char*, char*, Rule*)) R_GetCCallable("RxODE","new_string");
  return fun(g, s, e, r);
}

Elem *new_utf8_char(Grammar *g, char *s, char *e, Rule *r){
  static Elem *(*fun)(Grammar*, char*, char*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(Grammar*, char*, char*, Rule*)) R_GetCCallable("RxODE","new_utf8_char");
  return fun(g, s, e, r);
}

Elem *new_ident(char *s, char *e, Rule *r){
  static Elem *(*fun)(char*, char*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(char*, char*, Rule*)) R_GetCCallable("RxODE","new_ident");
  return fun(s, e, r);
}

void new_token(Grammar *g, char *s, char *e){
  static void (*fun)(Grammar*, char*, char*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, char*, char*)) R_GetCCallable("RxODE","new_token");
  return fun(g, s, e);
}

Elem *new_code(Grammar *g, char *s, char *e, Rule *r){
  static Elem *(*fun)(Grammar*, char*, char*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(Grammar*, char*, char*, Rule*)) R_GetCCallable("RxODE","new_code");
  return fun(g, s, e, r);
}

void add_global_code(Grammar *g, char *start, char *end, int line){
  static void (*fun)(Grammar*, char*, char*, int)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, char*, char*, int)) R_GetCCallable("RxODE","add_global_code");
  return fun(g, start, end, line);
}

Production *new_internal_production(Grammar *g, Production *p){
  static Production *(*fun)(Grammar*, Production*)=NULL;
  if (fun == NULL) fun = (Production* (*)(Grammar*, Production*)) R_GetCCallable("RxODE","new_internal_production");
  return fun(g, p);
}

Elem *dup_elem(Elem *e, Rule *r){
  static Elem *(*fun)(Elem*, Rule*)=NULL;
  if (fun == NULL) fun = (Elem* (*)(Elem*, Rule*)) R_GetCCallable("RxODE","dup_elem");
  return fun(e, r);
}

void add_declaration(Grammar *g, char *start, char *end, uint kind, uint line){
  static void (*fun)(Grammar*, char*, char*, uint, uint)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, char*, char*, uint, uint)) R_GetCCallable("RxODE","add_declaration");
  return fun(g, start, end, kind, line);
}

void add_pass(Grammar *g, char *start, char *end, uint kind, uint line){
  static void (*fun)(Grammar*, char*, char*, uint, uint)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, char*, char*, uint, uint)) R_GetCCallable("RxODE","add_pass");
  return fun(g, start, end, kind, line);
}

D_Pass *find_pass(Grammar *g, char *start, char *end){
  static D_Pass *(*fun)(Grammar*, char*, char*)=NULL;
  if (fun == NULL) fun = (D_Pass* (*)(Grammar*, char*, char*)) R_GetCCallable("RxODE","find_pass");
  return fun(g, start, end);
}

void   conditional_EBNF (Grammar *g ){
  static void   (*fun)(Grammar* )=NULL;
  if (fun == NULL) fun = (void   (*)(Grammar* )) R_GetCCallable("RxODE","conditional_EBNF ");
  return fun(g );
}

void   star_EBNF (Grammar *g ){
  static void   (*fun)(Grammar* )=NULL;
  if (fun == NULL) fun = (void   (*)(Grammar* )) R_GetCCallable("RxODE","star_EBNF ");
  return fun(g );
}

void   plus_EBNF (Grammar *g ){
  static void   (*fun)(Grammar* )=NULL;
  if (fun == NULL) fun = (void   (*)(Grammar* )) R_GetCCallable("RxODE","plus_EBNF ");
  return fun(g );
}

void rep_EBNF(Grammar *g, int minimum, int maximum){
  static void (*fun)(Grammar*, int, int)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*, int, int)) R_GetCCallable("RxODE","rep_EBNF");
  return fun(g, minimum, maximum);
}

void initialize_productions(Grammar *g){
  static void (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*)) R_GetCCallable("RxODE","initialize_productions");
  return fun(g);
}

int state_for_declaration(Grammar *g, int iproduction){
  static int (*fun)(Grammar*, int)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*, int)) R_GetCCallable("RxODE","state_for_declaration");
  return fun(g, iproduction);
}

void build_scanners(struct Grammar *g){
  static void (*fun)(struct Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(struct Grammar*)) R_GetCCallable("RxODE","build_scanners");
  return fun(g);
}

void build_LR_tables(Grammar *g){
  static void (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (void (*)(Grammar*)) R_GetCCallable("RxODE","build_LR_tables");
  return fun(g);
}

void sort_VecAction(VecAction *v){
  static void (*fun)(VecAction*)=NULL;
  if (fun == NULL) fun = (void (*)(VecAction*)) R_GetCCallable("RxODE","sort_VecAction");
  return fun(v);
}

uint elem_symbol(Grammar *g, Elem *e){
  static uint (*fun)(Grammar*, Elem*)=NULL;
  if (fun == NULL) fun = (uint (*)(Grammar*, Elem*)) R_GetCCallable("RxODE","elem_symbol");
  return fun(g, e);
}

State *goto_State(State *s, Elem *e){
  static State *(*fun)(State*, Elem*)=NULL;
  if (fun == NULL) fun = (State* (*)(State*, Elem*)) R_GetCCallable("RxODE","goto_State");
  return fun(s, e);
}

void free_Action(Action *a){
  static void (*fun)(Action*)=NULL;
  if (fun == NULL) fun = (void (*)(Action*)) R_GetCCallable("RxODE","free_Action");
  return fun(a);
}

void mkdparse(struct Grammar *g, char *grammar_pathname){
  static void (*fun)(struct Grammar*, char*)=NULL;
  if (fun == NULL) fun = (void (*)(struct Grammar*, char*)) R_GetCCallable("RxODE","mkdparse");
  return fun(g, grammar_pathname);
}

void mkdparse_from_string(struct Grammar *g, char *str){
  static void (*fun)(struct Grammar*, char*)=NULL;
  if (fun == NULL) fun = (void (*)(struct Grammar*, char*)) R_GetCCallable("RxODE","mkdparse_from_string");
  return fun(g, str);
}

D_ParseNode *ambiguity_count_fn(D_Parser *pp, int n, D_ParseNode **v){
  static D_ParseNode *(*fun)(D_Parser*, int, D_ParseNode**)=NULL;
  if (fun == NULL) fun = (D_ParseNode* (*)(D_Parser*, int, D_ParseNode**)) R_GetCCallable("RxODE","ambiguity_count_fn");
  return fun(pp, n, v);
}

BinaryTables *read_binary_tables(char *file_name, D_ReductionCode spec_code, D_ReductionCode final_code){
  static BinaryTables *(*fun)(char*, D_ReductionCode, D_ReductionCode)=NULL;
  if (fun == NULL) fun = (BinaryTables* (*)(char*, D_ReductionCode, D_ReductionCode)) R_GetCCallable("RxODE","read_binary_tables");
  return fun(file_name, spec_code, final_code);
}

BinaryTables *read_binary_tables_from_file(FILE *fp, D_ReductionCode spec_code, D_ReductionCode final_code){
  static BinaryTables *(*fun)(FILE*, D_ReductionCode, D_ReductionCode)=NULL;
  if (fun == NULL) fun = (BinaryTables* (*)(FILE*, D_ReductionCode, D_ReductionCode)) R_GetCCallable("RxODE","read_binary_tables_from_file");
  return fun(fp, spec_code, final_code);
}

BinaryTables *read_binary_tables_from_string(unsigned char *buf, D_ReductionCode spec_code, D_ReductionCode final_code){
  static BinaryTables *(*fun)(unsigned char*, D_ReductionCode, D_ReductionCode)=NULL;
  if (fun == NULL) fun = (BinaryTables* (*)(unsigned char*, D_ReductionCode, D_ReductionCode)) R_GetCCallable("RxODE","read_binary_tables_from_string");
  return fun(buf, spec_code, final_code);
}

void free_BinaryTables(BinaryTables * binary_tables){
  static void (*fun)(BinaryTables*)=NULL;
  if (fun == NULL) fun = (void (*)(BinaryTables*)) R_GetCCallable("RxODE","free_BinaryTables");
  return fun(binary_tables);
}

int scan_buffer(d_loc_t *loc, D_State *st, ShiftResult *result){
  static int (*fun)(d_loc_t*, D_State*, ShiftResult*)=NULL;
  if (fun == NULL) fun = (int (*)(d_loc_t*, D_State*, ShiftResult*)) R_GetCCallable("RxODE","scan_buffer");
  return fun(loc, st, result);
}

void vec_add_internal(void *v,  void *elem){
  static void (*fun)(void*, void*)=NULL;
  if (fun == NULL) fun = (void (*)(void*, void*)) R_GetCCallable("RxODE","vec_add_internal");
  return fun(v, elem);
}

int vec_eq(void *v, void *vv){
  static int (*fun)(void*, void*)=NULL;
  if (fun == NULL) fun = (int (*)(void*, void*)) R_GetCCallable("RxODE","vec_eq");
  return fun(v, vv);
}

int set_find(void *v, void *t){
  static int (*fun)(void*, void*)=NULL;
  if (fun == NULL) fun = (int (*)(void*, void*)) R_GetCCallable("RxODE","set_find");
  return fun(v, t);
}

int set_add(void *v, void *t){
  static int (*fun)(void*, void*)=NULL;
  if (fun == NULL) fun = (int (*)(void*, void*)) R_GetCCallable("RxODE","set_add");
  return fun(v, t);
}

int set_union(void *v, void *vv){
  static int (*fun)(void*, void*)=NULL;
  if (fun == NULL) fun = (int (*)(void*, void*)) R_GetCCallable("RxODE","set_union");
  return fun(v, vv);
}

void *set_add_fn(void *v, void *t, hash_fns_t *fns){
  static void *(*fun)(void*, void*, hash_fns_t*)=NULL;
  if (fun == NULL) fun = (void* (*)(void*, void*, hash_fns_t*)) R_GetCCallable("RxODE","set_add_fn");
  return fun(v, t, fns);
}

void set_union_fn(void *v, void *vv, hash_fns_t *fns){
  static void (*fun)(void*, void*, hash_fns_t*)=NULL;
  if (fun == NULL) fun = (void (*)(void*, void*, hash_fns_t*)) R_GetCCallable("RxODE","set_union_fn");
  return fun(v, vv, fns);
}

void set_to_vec(void *av){
  static void (*fun)(void*)=NULL;
  if (fun == NULL) fun = (void (*)(void*)) R_GetCCallable("RxODE","set_to_vec");
  return fun(av);
}

int buf_read(const char *pathname, char **buf, int *len){
  static int (*fun)(const char*, char**, int*)=NULL;
  if (fun == NULL) fun = (int (*)(const char*, char**, int*)) R_GetCCallable("RxODE","buf_read");
  return fun(pathname, buf, len);
}

char *sbuf_read(const char *pathname){
  static char *(*fun)(const char*)=NULL;
  if (fun == NULL) fun = (char* (*)(const char*)) R_GetCCallable("RxODE","sbuf_read");
  return fun(pathname);
}

char *dup_str(const char *str, const char *end){
  static char *(*fun)(const char*, const char*)=NULL;
  if (fun == NULL) fun = (char* (*)(const char*, const char*)) R_GetCCallable("RxODE","dup_str");
  return fun(str, end);
}

uint strhashl(const char *s, int len){
  static uint (*fun)(const char*, int)=NULL;
  if (fun == NULL) fun = (uint (*)(const char*, int)) R_GetCCallable("RxODE","strhashl");
  return fun(s, len);
}

void d_free(void *x){
  static void (*fun)(void*)=NULL;
  if (fun == NULL) fun = (void (*)(void*)) R_GetCCallable("RxODE","d_free");
  return fun(x);
}

void int_list_diff(int *a, int *b, int *c){
  static void (*fun)(int*, int*, int*)=NULL;
  if (fun == NULL) fun = (void (*)(int*, int*, int*)) R_GetCCallable("RxODE","int_list_diff");
  return fun(a, b, c);
}

void int_list_intersect(int *a, int *b, int *c){
  static void (*fun)(int*, int*, int*)=NULL;
  if (fun == NULL) fun = (void (*)(int*, int*, int*)) R_GetCCallable("RxODE","int_list_intersect");
  return fun(a, b, c);
}

int *int_list_dup(int *aa){
  static int *(*fun)(int*)=NULL;
  if (fun == NULL) fun = (int* (*)(int*)) R_GetCCallable("RxODE","int_list_dup");
  return fun(aa);
}

char *escape_string(char *s){
  static char *(*fun)(char*)=NULL;
  if (fun == NULL) fun = (char* (*)(char*)) R_GetCCallable("RxODE","escape_string");
  return fun(s);
}

char *escape_string_single_quote(char *s){
  static char *(*fun)(char*)=NULL;
  if (fun == NULL) fun = (char* (*)(char*)) R_GetCCallable("RxODE","escape_string_single_quote");
  return fun(s);
}

int write_c_tables(Grammar *g){
  static int (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*)) R_GetCCallable("RxODE","write_c_tables");
  return fun(g);
}

int write_binary_tables(Grammar *g){
  static int (*fun)(Grammar*)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*)) R_GetCCallable("RxODE","write_binary_tables");
  return fun(g);
}

int write_binary_tables_to_file(Grammar *g, FILE *fp){
  static int (*fun)(Grammar*, FILE*)=NULL;
  if (fun == NULL) fun = (int (*)(Grammar*, FILE*)) R_GetCCallable("RxODE","write_binary_tables_to_file");
  return fun(g, fp);
}

void set_d_use_r_headers(int x){
  static void (*fun)(int)=NULL;
  if (fun == NULL) fun = (void (*)(int)) R_GetCCallable("RxODE","set_d_use_r_headers");
  return fun(x);
}

void set_d_rdebug_grammar_level(int x){
  static void (*fun)(int)=NULL;
  if (fun == NULL) fun = (void (*)(int)) R_GetCCallable("RxODE","set_d_rdebug_grammar_level");
  return fun(x);
}

void set_d_use_file_name(int x){
  static void (*fun)(int)=NULL;
  if (fun == NULL) fun = (void (*)(int)) R_GetCCallable("RxODE","set_d_use_file_name");
  return fun(x);
}

void set_d_verbose_level(int x){
  static void (*fun)(int)=NULL;
  if (fun == NULL) fun = (void (*)(int)) R_GetCCallable("RxODE","set_d_verbose_level");
  return fun(x);
}

void set_d_debug_level(int x){
  static void (*fun)(int)=NULL;
  if (fun == NULL) fun = (void (*)(int)) R_GetCCallable("RxODE","set_d_debug_level");
  return fun(x);
}

int get_d_use_r_headers(){
  static int (*fun)()=NULL;
  if (fun == NULL) fun = (int (*)()) R_GetCCallable("RxODE","get_d_use_r_headers");
  return fun();
}

int get_d_rdebug_grammar_level(){
  static int (*fun)()=NULL;
  if (fun == NULL) fun = (int (*)()) R_GetCCallable("RxODE","get_d_rdebug_grammar_level");
  return fun();
}

int get_d_use_file_name(){
  static int (*fun)()=NULL;
  if (fun == NULL) fun = (int (*)()) R_GetCCallable("RxODE","get_d_use_file_name");
  return fun();
}

int get_d_verbose_level(){
  static int (*fun)()=NULL;
  if (fun == NULL) fun = (int (*)()) R_GetCCallable("RxODE","get_d_verbose_level");
  return fun();
}

int get_d_debug_level(){
  static int (*fun)()=NULL;
  if (fun == NULL) fun = (int (*)()) R_GetCCallable("RxODE","get_d_debug_level");
  return fun();
}

void set_d_file_name(char *x){
  static void (*fun)(char*)=NULL;
  if (fun == NULL) fun = (void (*)(char*)) R_GetCCallable("RxODE","set_d_file_name");
  return fun(x);
}

#if defined(__cplusplus)
}
#endif
#endif
