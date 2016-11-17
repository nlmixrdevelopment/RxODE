/*
Register C callables to R.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Rdynload.h>
#include "gramgram.h"
#include "d.h"
#include "mkdparse.h"
#include "dparse.h"

extern int d_use_r_headers;
void set_d_use_r_headers(int x){
  d_use_r_headers = x;
}

extern int d_rdebug_grammar_level;
void set_d_rdebug_grammar_level(int x){
  d_rdebug_grammar_level = x;
}

extern int d_use_file_name;
void set_d_use_file_name(int x){
  d_use_file_name = x;
}

extern int d_verbose_level;
void set_d_verbose_level(int x){
  d_verbose_level = x;
}

extern int d_debug_level;
void set_d_debug_level(int x){
  d_debug_level = x;
}


extern char * d_file_name;
void set_d_file_name(char *x){
  d_file_name=x;
}


void R_ini_RxODE(DllInfo *info);

D_Parser *new_D_Parser(struct D_ParserTables *t, int sizeof_ParseNode_User);
void free_D_Parser(D_Parser *p); 
D_ParseNode *dparse(D_Parser *p, char *buf, int buf_len);
void free_D_ParseNode(D_Parser *p, D_ParseNode *pn);
void free_D_ParseTreeBelow(D_Parser *p, D_ParseNode *pn);
int d_get_number_of_children(D_ParseNode *pn);
D_ParseNode *d_get_child(D_ParseNode *pn, int child);
D_ParseNode *d_find_in_tree(D_ParseNode *pn, int symbol);
char *d_ws_before(D_Parser *p, D_ParseNode *pn); 
char *d_ws_after(D_Parser *p, D_ParseNode *pn); 
void d_pass(D_Parser *p, D_ParseNode *pn, int pass_number);
int resolve_amb_greedy(D_Parser *dp, int n, D_ParseNode **v);
char *d_dup_pathname_str(const char *str);
void parse_whitespace(struct D_Parser *p, d_loc_t *loc, void **p_globals);
D_Scope *new_D_Scope(D_Scope *parent);
D_Scope *enter_D_Scope(D_Scope *current, D_Scope *scope);
D_Scope *commit_D_Scope(D_Scope *scope);
D_Scope *equiv_D_Scope(D_Scope *scope);
D_Scope *global_D_Scope(D_Scope *scope);
D_Scope *scope_D_Scope(D_Scope *current, D_Scope *scope);
void free_D_Scope(D_Scope *st, int force);
D_Sym *new_D_Sym(D_Scope *st, char *name, char *end, int sizeof_D_Sym);
void free_D_Sym(D_Sym *sym);
D_Sym *find_D_Sym(D_Scope *st, char *name, char *end);
D_Sym *find_global_D_Sym(D_Scope *st, char *name, char *end);
D_Sym *update_D_Sym(D_Sym *sym, D_Scope **st, int sizeof_D_Sym);
D_Sym *update_additional_D_Sym(D_Scope *st, D_Sym *sym, int sizeof_D_Sym);
D_Sym *current_D_Sym(D_Scope *st, D_Sym *sym);
D_Sym *find_D_Sym_in_Scope(D_Scope *st, D_Scope *cur, char *name, char *end);
D_Sym *next_D_Sym_in_Scope(D_Scope **st, D_Sym **sym);
void print_scope(D_Scope *st);
Grammar *new_D_Grammar(char *pathname);
void free_D_Grammar(Grammar *g);
int build_grammar(Grammar *g);
int parse_grammar(Grammar *g, char *pathname, char *str);
void print_grammar(Grammar *g);
void print_rdebug_grammar(Grammar *g, char *pathname);
void print_states(Grammar *g);
void print_rule(Rule *r);
void print_term(Term *t);
Production *lookup_production(Grammar *g, char *name, int len);
Rule *new_rule(Grammar *g, Production *p);
Elem *new_elem_nterm(Production *p, Rule *r);
void new_declaration(Grammar *g, Elem *e, uint kind);
Production *new_production(Grammar *g, char *name);
Elem *new_string(Grammar *g, char *s, char *e, Rule *r);
Elem *new_utf8_char(Grammar *g, char *s, char *e, Rule *r);
Elem *new_ident(char *s, char *e, Rule *r);
void new_token(Grammar *g, char *s, char *e);
Elem *new_code(Grammar *g, char *s, char *e, Rule *r);
void add_global_code(Grammar *g, char *start, char *end, int line);
Production *new_internal_production(Grammar *g, Production *p);
Elem * dup_elem(Elem *e, Rule *r);
void add_declaration(Grammar *g, char *start, char *end, uint kind, uint line);
void add_pass(Grammar *g, char *start, char *end, uint kind, uint line);
D_Pass *find_pass(Grammar *g, char *start, char *end);
void conditional_EBNF(Grammar *g); 
void star_EBNF(Grammar *g); 
void plus_EBNF(Grammar *g); 
void rep_EBNF(Grammar *g, int minimum, int maximum);
void initialize_productions(Grammar *g);
int state_for_declaration(Grammar *g, int iproduction);
void build_scanners(struct Grammar *g);
void build_LR_tables(Grammar *g);
void sort_VecAction(VecAction *v);
uint elem_symbol(Grammar *g, Elem *e);
State *goto_State(State *s, Elem *e);
void free_Action(Action *a);
void mkdparse(struct Grammar *g, char *grammar_pathname);
void mkdparse_from_string(struct Grammar *g, char *str);
D_ParseNode * ambiguity_count_fn(D_Parser *pp, int n, D_ParseNode **v);
BinaryTables * read_binary_tables(char *file_name, D_ReductionCode spec_code, D_ReductionCode final_code);
BinaryTables * read_binary_tables_from_file(FILE *fp, D_ReductionCode spec_code, D_ReductionCode final_code);
BinaryTables * read_binary_tables_from_string(unsigned char *buf, D_ReductionCode spec_code, D_ReductionCode final_code);
void free_BinaryTables(BinaryTables * binary_tables);
int scan_buffer(d_loc_t *loc, D_State *st, ShiftResult *result);
void vec_add_internal(void *v,  void *elem);
int vec_eq(void *v, void *vv);
int set_find(void *v, void *t);
int set_add(void *v, void *t);
int set_union(void *v, void *vv);
void *set_add_fn(void *v, void *t, hash_fns_t *fns);
void set_union_fn(void *v, void *vv, hash_fns_t *fns);
void set_to_vec(void *av);
void * stack_push_internal(AbstractStack*, void*);
int buf_read(const char *pathname, char **buf, int *len);
char *sbuf_read(const char *pathname);
char *dup_str(const char *str, const char *end);
uint strhashl(const char *s, int len);
void d_free(void *);
void int_list_diff(int *a, int *b, int *c);
void int_list_intersect(int *a, int *b, int *c);
int *int_list_dup(int *aa);
char *escape_string(char *s);
char *escape_string_single_quote(char *s);
int write_c_tables(Grammar *g);
int write_binary_tables(Grammar *g);
int write_binary_tables_to_file(Grammar *g, FILE *fp);
void set_d_use_r_headers(int x);
void set_d_rdebug_grammar_level(int x);
void set_d_use_file_name(int x);
void set_d_verbose_level(int x);
void set_d_debug_level(int x);
void set_d_file_name(char *x);

void R_init_RxODE(DllInfo *info){
  R_ini_RxODE(info);
  R_RegisterCCallable("RxODE","set_d_file_name",(DL_FUNC) set_d_file_name);
  R_RegisterCCallable("RxODE","set_d_debug_level",(DL_FUNC) set_d_debug_level);
  R_RegisterCCallable("RxODE","set_d_verbose_level",(DL_FUNC) set_d_verbose_level);
  R_RegisterCCallable("RxODE","set_d_use_file_name",(DL_FUNC) set_d_use_file_name);
  R_RegisterCCallable("RxODE","set_d_rdebug_grammar_level",(DL_FUNC) set_d_rdebug_grammar_level);
  R_RegisterCCallable("RxODE","set_d_use_r_headers",(DL_FUNC) set_d_use_r_headers);
  R_RegisterCCallable("RxODE","write_binary_tables_to_file",(DL_FUNC) write_binary_tables_to_file);
  R_RegisterCCallable("RxODE","write_binary_tables",(DL_FUNC) write_binary_tables);
  R_RegisterCCallable("RxODE","write_c_tables",(DL_FUNC) write_c_tables);
  R_RegisterCCallable("RxODE","escape_string_single_quote",(DL_FUNC) escape_string_single_quote);
  R_RegisterCCallable("RxODE","escape_string",(DL_FUNC) escape_string);
  R_RegisterCCallable("RxODE","int_list_dup",(DL_FUNC) int_list_dup);
  R_RegisterCCallable("RxODE","int_list_intersect",(DL_FUNC) int_list_intersect);
  R_RegisterCCallable("RxODE","int_list_diff",(DL_FUNC) int_list_diff);
  R_RegisterCCallable("RxODE","d_free",(DL_FUNC) d_free);
  R_RegisterCCallable("RxODE","strhashl",(DL_FUNC) strhashl);
  R_RegisterCCallable("RxODE","dup_str",(DL_FUNC) dup_str);
  R_RegisterCCallable("RxODE","sbuf_read",(DL_FUNC) sbuf_read);
  R_RegisterCCallable("RxODE","buf_read",(DL_FUNC) buf_read);
  R_RegisterCCallable("RxODE","stack_push_internal",(DL_FUNC) stack_push_internal);
  R_RegisterCCallable("RxODE","set_to_vec",(DL_FUNC) set_to_vec);
  R_RegisterCCallable("RxODE","set_union_fn",(DL_FUNC) set_union_fn);
  R_RegisterCCallable("RxODE","set_add_fn",(DL_FUNC) set_add_fn);
  R_RegisterCCallable("RxODE","set_union",(DL_FUNC) set_union);
  R_RegisterCCallable("RxODE","set_add",(DL_FUNC) set_add);
  R_RegisterCCallable("RxODE","set_find",(DL_FUNC) set_find);
  R_RegisterCCallable("RxODE","vec_eq",(DL_FUNC) vec_eq);
  R_RegisterCCallable("RxODE","vec_add_internal",(DL_FUNC) vec_add_internal);
  R_RegisterCCallable("RxODE","scan_buffer",(DL_FUNC) scan_buffer);
  R_RegisterCCallable("RxODE","free_BinaryTables",(DL_FUNC) free_BinaryTables);
  R_RegisterCCallable("RxODE","read_binary_tables_from_string",(DL_FUNC) read_binary_tables_from_string);
  R_RegisterCCallable("RxODE","read_binary_tables_from_file",(DL_FUNC) read_binary_tables_from_file);
  R_RegisterCCallable("RxODE","read_binary_tables",(DL_FUNC) read_binary_tables);
  R_RegisterCCallable("RxODE","ambiguity_count_fn",(DL_FUNC) ambiguity_count_fn);
  R_RegisterCCallable("RxODE","mkdparse_from_string",(DL_FUNC) mkdparse_from_string);
  R_RegisterCCallable("RxODE","mkdparse",(DL_FUNC) mkdparse);
  R_RegisterCCallable("RxODE","free_Action",(DL_FUNC) free_Action);
  R_RegisterCCallable("RxODE","goto_State",(DL_FUNC) goto_State);
  R_RegisterCCallable("RxODE","elem_symbol",(DL_FUNC) elem_symbol);
  R_RegisterCCallable("RxODE","sort_VecAction",(DL_FUNC) sort_VecAction);
  R_RegisterCCallable("RxODE","build_LR_tables",(DL_FUNC) build_LR_tables);
  R_RegisterCCallable("RxODE","build_scanners",(DL_FUNC) build_scanners);
  R_RegisterCCallable("RxODE","state_for_declaration",(DL_FUNC) state_for_declaration);
  R_RegisterCCallable("RxODE","initialize_productions",(DL_FUNC) initialize_productions);
  R_RegisterCCallable("RxODE","rep_EBNF",(DL_FUNC) rep_EBNF);
  R_RegisterCCallable("RxODE","plus_EBNF ",(DL_FUNC) plus_EBNF );
  R_RegisterCCallable("RxODE","star_EBNF ",(DL_FUNC) star_EBNF );
  R_RegisterCCallable("RxODE","conditional_EBNF ",(DL_FUNC) conditional_EBNF );
  R_RegisterCCallable("RxODE","find_pass",(DL_FUNC) find_pass);
  R_RegisterCCallable("RxODE","add_pass",(DL_FUNC) add_pass);
  R_RegisterCCallable("RxODE","add_declaration",(DL_FUNC) add_declaration);
  R_RegisterCCallable("RxODE","dup_elem",(DL_FUNC) dup_elem);
  R_RegisterCCallable("RxODE","new_internal_production",(DL_FUNC) new_internal_production);
  R_RegisterCCallable("RxODE","add_global_code",(DL_FUNC) add_global_code);
  R_RegisterCCallable("RxODE","new_code",(DL_FUNC) new_code);
  R_RegisterCCallable("RxODE","new_token",(DL_FUNC) new_token);
  R_RegisterCCallable("RxODE","new_ident",(DL_FUNC) new_ident);
  R_RegisterCCallable("RxODE","new_utf8_char",(DL_FUNC) new_utf8_char);
  R_RegisterCCallable("RxODE","new_string",(DL_FUNC) new_string);
  R_RegisterCCallable("RxODE","new_production",(DL_FUNC) new_production);
  R_RegisterCCallable("RxODE","new_declaration",(DL_FUNC) new_declaration);
  R_RegisterCCallable("RxODE","new_elem_nterm",(DL_FUNC) new_elem_nterm);
  R_RegisterCCallable("RxODE","new_rule",(DL_FUNC) new_rule);
  R_RegisterCCallable("RxODE","lookup_production",(DL_FUNC) lookup_production);
  R_RegisterCCallable("RxODE","print_term",(DL_FUNC) print_term);
  R_RegisterCCallable("RxODE","print_rule",(DL_FUNC) print_rule);
  R_RegisterCCallable("RxODE","print_states",(DL_FUNC) print_states);
  R_RegisterCCallable("RxODE","print_rdebug_grammar",(DL_FUNC) print_rdebug_grammar);
  R_RegisterCCallable("RxODE","print_grammar",(DL_FUNC) print_grammar);
  R_RegisterCCallable("RxODE","parse_grammar",(DL_FUNC) parse_grammar);
  R_RegisterCCallable("RxODE","build_grammar",(DL_FUNC) build_grammar);
  R_RegisterCCallable("RxODE","free_D_Grammar",(DL_FUNC) free_D_Grammar);
  R_RegisterCCallable("RxODE","new_D_Grammar",(DL_FUNC) new_D_Grammar);
  R_RegisterCCallable("RxODE","print_scope",(DL_FUNC) print_scope);
  R_RegisterCCallable("RxODE","next_D_Sym_in_Scope",(DL_FUNC) next_D_Sym_in_Scope);
  R_RegisterCCallable("RxODE","find_D_Sym_in_Scope",(DL_FUNC) find_D_Sym_in_Scope);
  R_RegisterCCallable("RxODE","current_D_Sym",(DL_FUNC) current_D_Sym);
  R_RegisterCCallable("RxODE","update_additional_D_Sym",(DL_FUNC) update_additional_D_Sym);
  R_RegisterCCallable("RxODE","update_D_Sym",(DL_FUNC) update_D_Sym);
  R_RegisterCCallable("RxODE","find_global_D_Sym",(DL_FUNC) find_global_D_Sym);
  R_RegisterCCallable("RxODE","find_D_Sym",(DL_FUNC) find_D_Sym);
  R_RegisterCCallable("RxODE","free_D_Sym",(DL_FUNC) free_D_Sym);
  R_RegisterCCallable("RxODE","new_D_Sym",(DL_FUNC) new_D_Sym);
  R_RegisterCCallable("RxODE","free_D_Scope",(DL_FUNC) free_D_Scope);
  R_RegisterCCallable("RxODE","scope_D_Scope",(DL_FUNC) scope_D_Scope);
  R_RegisterCCallable("RxODE","global_D_Scope",(DL_FUNC) global_D_Scope);
  R_RegisterCCallable("RxODE","equiv_D_Scope",(DL_FUNC) equiv_D_Scope);
  R_RegisterCCallable("RxODE","commit_D_Scope",(DL_FUNC) commit_D_Scope);
  R_RegisterCCallable("RxODE","enter_D_Scope",(DL_FUNC) enter_D_Scope);
  R_RegisterCCallable("RxODE","new_D_Scope",(DL_FUNC) new_D_Scope);
  R_RegisterCCallable("RxODE","parse_whitespace",(DL_FUNC) parse_whitespace);
  R_RegisterCCallable("RxODE","d_dup_pathname_str",(DL_FUNC) d_dup_pathname_str);
  R_RegisterCCallable("RxODE","resolve_amb_greedy",(DL_FUNC) resolve_amb_greedy);
  R_RegisterCCallable("RxODE","d_pass",(DL_FUNC) d_pass);
  R_RegisterCCallable("RxODE","d_ws_after ",(DL_FUNC) d_ws_after );
  R_RegisterCCallable("RxODE","d_ws_before ",(DL_FUNC) d_ws_before );
  R_RegisterCCallable("RxODE","d_find_in_tree",(DL_FUNC) d_find_in_tree);
  R_RegisterCCallable("RxODE","d_get_child",(DL_FUNC) d_get_child);
  R_RegisterCCallable("RxODE","d_get_number_of_children",(DL_FUNC) d_get_number_of_children);
  R_RegisterCCallable("RxODE","free_D_ParseTreeBelow",(DL_FUNC) free_D_ParseTreeBelow);
  R_RegisterCCallable("RxODE","free_D_ParseNode",(DL_FUNC) free_D_ParseNode);
  R_RegisterCCallable("RxODE","dparse",(DL_FUNC) dparse);
  R_RegisterCCallable("RxODE","free_D_Parser ",(DL_FUNC) free_D_Parser );
  R_RegisterCCallable("RxODE","new_D_Parser",(DL_FUNC) new_D_Parser);

}

