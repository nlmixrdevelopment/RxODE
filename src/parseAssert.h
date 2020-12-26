////////////////////////////////////////////////////////////////////////////////
// assertions
static inline int assertNoRAssign(nodeInfo ni, char *name, D_ParseNode *pn, int i){
  if (!rx_syntax_assign  &&
      ((i == 4 && nodeHas(derivative)) ||
       (i == 6 && nodeHas(dfdy)))) {
    D_ParseNode *xpn = d_get_child(pn,i);
    char *v = (char*)rc_dup_str(xpn->start_loc.s, xpn->end);
    if (!strcmp("<-",v)){
      updateSyntaxCol();
      trans_syntax_error_report_fn(NOASSIGN);
    }
    /* Free(v); */
    /* continue; */
    return 1;
  }
  return 0;
}

static inline void assertEndSemicolon(nodeInfo ni, char *name, int i, D_ParseNode *xpn) {
  if (rx_syntax_require_semicolon && nodeHas(end_statement) && i == 0){
    if (xpn->start_loc.s ==  xpn->end){
      updateSyntaxCol();
      trans_syntax_error_report_fn(NEEDSEMI);
    }
  }
}
