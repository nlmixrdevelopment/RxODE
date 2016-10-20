##' Make a dparser c file based on a grammer 
##'
##' Uses a grammer file to create a c file for parsing.
##'
##' This is used internally to create the parser in trans
##'
##'  mkdparser is a scannerless GLR parser generator based on the
##'  Tomita algorithm. It is self-hosted and very easy to
##'  use. Grammars are written in a natural style of EBNF and regular
##'  expressions and support both speculative and final actions.
##' 
##' @title mkdparse dparser grammer c
##' 
##' @param file File name of grammer to parse.
##' 
##' @param outputFile Output file name.  Can be a directory.  If so,
##'     the file name is determined by the input file.
##' 
##' @param set_op_priority_from_rule Toggle setting of operator
##'     priority from rules.  Setting of operator priorities on
##'     operator tokens can increase the size of the tables but can
##'     permit unnecessary parse stacks to be pruned earlier. (FALSE
##'     by default)
##' 
##' @param right_recursive_BNF Toggle use of right recursion for EBNF
##'     productions.  Do not change this unless you really know what
##'     you are doing. (FALSE by default)
##' 
##' @param states_for_whitespace Toggle computing whitespace states.
##'     If 'whitespace' is defined in the grammar, then use it as a
##'     subparser to consume whitespace. (TRUE by default)
##' 
##' @param states_for_all_nterms Toggle computing states for all
##'     non-terminals.  Ensures that there is a unique state for each
##'     non-terminal so that a subparsers can be invoked for that
##'     non-terminal. (FALSE by default)
##' 
##' @param tokenizer Toggle building of a tokenizer for START.  When
##'     TRUE, instead of generating a unique scanner for each state
##'     (i.e. a 'scannerless' parser), generate a single scanner
##'     (tokenizer) for the entire grammar.  This provides an easy way
##'     to build grammars for languages which assume a tokenizer
##'     (e.g. ANSI C). (FALSE by default)
##' 
##' @param token_type Token type "#define" or "enum"
##' 
##' @param longest_match Toggle longest match lexical ambiguity
##'     resolution.  When TRUE the scanner only recognizing the
##'     longest matching tokens in a given state. This provides an
##'     easy way to build grammars for languages which use longest
##'     match lexical ambiguity resolution (e.g. ANSI-C, C++). (FALSE
##'     by default)
##' 
##' @param grammar_ident Tag for grammar data structures so that
##'     multiple sets of tables can be included in one
##'     file/application. (defaults to 'gram')
##' 
##' @param ident_from_filename Build the grammer tag from the
##'     file-name.
##' 
##' @param scanner_blocks Number of blocks to which scanner tables are
##'     broken up into. Larger numbers permit more sharing with more
##'     overhead.  4 seems to be optimal for most grammars. (defaults
##'     to 4) files.
##' 
##' @param write_line_directives Toggle writing of line numbers.  Used
##'     to debug the parsing table generator itself. (TRUE by default)
##' 
##' @param write_header Write header, FALSE : no, TRUE : yes,
##'     "IfEmpty" : only if not empty.
##' 
##' @param rdebug Replace all actions in the grammar with actions
##'     printing productions, 1 : during the speculative parsing
##'     process (<-), 2 : when reduction is part of any legal final
##'     parse (<=), 3 : both, 4 : remove all actions from the grammar.
##'     Print the changed grammar to console.  Useful for debugging or
##'     prototyping new, experimental grammars.
##' 
##' @param verbose Increase verbosity.
##' 
##' @param write_extension Set the extension of the generated code
##'     file.  For C++ programs (for example) the extension can be set
##'     to .cpp with the option \code{write_extension="cpp"}.
##'     (\code{write_extension="c"} by default)
##' 
##' @return Nothing. Outputs files instead.
##' 
##' @author Matthew L. Fidler for R interface, John Plevyak for
##'     dparser
##' 
##' @keywords internal
##' @export
mkdparse <- function(file,outputFile,
                     set_op_priority_from_rule = FALSE,
                     right_recursive_BNF = FALSE,
                     states_for_whitespace = TRUE,
                     states_for_all_nterms = FALSE,
                     tokenizer = FALSE,
                     token_type = c("#define","enum"),
                     longest_match = FALSE,
                     grammar_ident ="gram",
                     ident_from_filename = FALSE,
                     scanner_blocks = 4,
                     write_line_directives = TRUE,
                     write_header = c("IfEmpty",TRUE,FALSE),
                     rdebug = FALSE,
                     verbose = TRUE,
                     write_extension="c"
                     ){
    if (missing(write_header) || write_header == "IfEmpty"){
        write_header <- -1;
    }
    if (ident_from_filename){
        ## Put ident from the filename
        grammar_ident = gsub("[.][^.]$","",basename(file));
    }
    if (missing(outputFile)){
        outputFile <- file.path(dirname(file),paste0(basename(file),".d_parser.",write_extension));
    } else if (dir.exists(outputFile)){
        outputFile <- file.path(outputFile,paste0(basename(file),".d_parser.",write_extension));
    }
    if (missing(token_type)){
        token_type <- 0;
    } else if (token_type == "enum"){
        token_type <- 1;
    } else {
        token_type <- 0;
    }
    .Call(cDparser,
          file,
          outputFile,
          as.integer(set_op_priority_from_rule),
          as.integer(right_recursive_BNF),
          as.integer(states_for_whitespace),
          as.integer(states_for_all_nterms),
          as.integer(tokenizer),
          as.integer(longest_match),
          grammar_ident,
          as.integer(scanner_blocks),
          as.integer(write_line_directives),
          as.integer(rdebug),
          as.integer(verbose),
          write_extension,
          as.integer(write_header),
          as.integer(token_type));
    return(invisible());
}

updateParser <- function(){
    mkdparse(devtools::package_file("inst/tran.g"),
             devtools::package_file("src/"));
}
