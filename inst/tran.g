//loop

statement_list : (statement)+ ;

statement 
  : ini end_statement 
  | ini0 end_statement
  | ini0f end_statement
  | assignment end_statement
  | derivative end_statement
  | jac end_statement
  | dfdy end_statement
  | compound_statement
  | selection_statement
  | printf_statement end_statement
  | print_command end_statement
  | end_statement ;


compound_statement : '{' statement_list? '}' ;

selection_statement
  : 'if' '(' logical_or_expression ')' statement ('else' statement)?;

printf_statement
  : printf_command '(' string (',' additive_expression )* ')';

printf_command
  : 'printf'      | 'Rprintf'      | 'print'      |
    'jac_printf'  | 'jac_Rprintf'  | 'jac_print'  |
    'ode_printf'  | 'ode_Rprintf'  | 'ode_print'  |
    'jac0_printf' | 'jac0_Rprintf' | 'jac0_print' |
    'ode_printf'  | 'ode_Rprintf'  | 'ode_print'  |
    'ode0_printf' | 'ode0_Rprintf' | 'ode0_print' |
    'lhs_printf'  | 'lhs_Rprintf'  | 'lhs_print'  ;

print_command
  : 'print' | 'ode_print' | 'jac_print' | 'lhs_print';

ini0       : identifier_r ('(0)' | '{0}' | '[0]') ('=' | '<-' ) ini_const;

ini0f     : identifier_r ('(0)' | '{0}' | '[0]') ('=' | '<-' ) additive_expression;

ini        : identifier_r ('=' | '<-' ) ini_const;

derivative : 'd/dt' '(' identifier_r_no_output ')' ('=' | '<-' | '~')  ('+' | '-' | ) additive_expression;
der_rhs    : 'd/dt' '(' identifier_r_no_output ')';
jac        : jac_command '(' identifier_r_no_output ',' (theta0_noout | theta_noout | eta_noout | identifier_r_no_output) ')' ('=' | '<-' ) additive_expression;
jac_rhs    : jac_command '(' identifier_r_no_output ',' (theta0_noout | theta_noout | eta_noout | identifier_r_no_output) ')';

inf_cmd    : ('rxRate');
inf_rhs    : inf_cmd '(' identifier_r_no_output ')';

// transit(n,mtt) -> transit3(t,n,mtt)
transit2   : 'transit' '(' trans_const ',' trans_const ')';

// transit(n,mtt, bio) -> transit4(t,n,mtt,bio)
transit3   : 'transit' '(' trans_const ',' trans_const ',' trans_const ')';

dfdy        : 'df' '(' identifier_r_no_output ')/dy(' (theta0_noout | theta_noout | eta_noout | identifier_r_no_output) ')' ('=' | '<-' ) additive_expression;
dfdy_rhs    : 'df' '(' identifier_r_no_output ')/dy(' (theta0_noout | theta_noout | eta_noout | identifier_r_no_output) ')';

jac_command : 'jac' | 'df/dy';

end_statement : (';')* ;

assignment : identifier_r ('=' | '<-' | '~' ) additive_expression;

logical_or_expression : logical_and_expression 
    (('||' | '|')  logical_and_expression)* ;

logical_and_expression : equality_expression0 
    (('&&' | '&') equality_expression0)* ;

equality_expression0 : equality_expression |
    '(' equality_expression ')' |
    '!' '(' equality_expression ')';

equality_expression : relational_expression 
    (('!=' | '~=' | '<>' | '/=' | '==') relational_expression)* ;

relational_expression : additive_expression
    (('<' | '>' | '<=' | '>=') additive_expression)* ;

additive_expression : multiplicative_expression
    (('+' | '-') multiplicative_expression)* ;

multiplicative_expression : unary_expression 
    (mult_part)* ;

mult_part : ('*' | '/') unary_expression ;

unary_expression : ('+' | '-')? (theta0 | theta | eta | primary_expression | power_expression );

exponent_expression : ('+' | '-')? (theta0 | theta | eta | primary_expression );

power_expression : primary_expression power_operator exponent_expression;

power_operator   : ('^' | '**');

primary_expression 
  : constant
  | identifier_r
  | theta0
  | theta
  | eta
  | der_rhs
  | inf_rhs
  | jac_rhs
  | dfdy_rhs
  | transit2
  | transit3
  | function
  | factorial 
  | factorial_exp 
  | lfactorial 
  | lfactorial_exp
  | '(' additive_expression ')'
  ;

factorial_exp: '(' additive_expression ')' '!';
factorial: fact_const '!';
lfactorial: "log *[(]" fact_const "! *[)]";
lfactorial_exp: "log *[(] *[(] *" additive_expression "[)] *! *[)]";

function : identifier '(' (additive_expression)* (',' additive_expression)* ')' ;

ini_const : '-'? constant;
fact_const: identifier_r_no_output | '-'? (decimalint | float1 | float2) ;
trans_const: identifier_r | '-'? constant;

constant : decimalint | float1 | float2;

identifier_r: identifier_r_1 | identifier_r_2 | 'transit';

identifier_r_no_output: identifier_r_no_output_1 | identifier_r_no_output_2 | 'transit';

theta: ('THETA' | 'theta') '[' decimalint ']';
eta: ('ETA' | 'eta') '[' decimalint ']';
theta0: ('THETA' | 'theta' | 'ETA' | 'eta');

theta_noout: ('THETA' | 'theta') '[' decimalint ']';
eta_noout: ('ETA' | 'eta') '[' decimalint ']';
theta0_noout: ('THETA' | 'theta' | 'ETA' | 'eta');


decimalint: "0|([1-9][0-9]*)" $term -1;
string: "\"([^\"\\]|\\[^])*\"";
float1: "([0-9]+.[0-9]*|[0-9]*.[0-9]+)([eE][\-\+]?[0-9]+)?" $term -2;
float2: "[0-9]+[eE][\-\+]?[0-9]+" $term -3;
identifier_r_1: "[a-zA-Z][a-zA-Z0-9_.]*" $term -4;
identifier_r_no_output_1: "[a-zA-Z][a-zA-Z0-9_.]*" $term -4;
identifier_r_2: "[.]+[a-zA-Z_][a-zA-Z0-9_.]*" $term -4;
identifier_r_no_output_2: "[.]+[a-zA-Z_][a-zA-Z0-9_.]*" $term -4;
identifier: "[a-zA-Z][a-zA-Z0-9_]*" $term -4;
whitespace: ( "[ \t\r\n]+" | singleLineComment )*;
singleLineComment: '#' "[^\n]*" '\n';

