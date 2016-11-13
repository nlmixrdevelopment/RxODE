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
##' @param use_r_header when TRUE, add R headers and swap printf for
##'     Rprintf. By default this is FALSE.
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
                     write_extension="c",
                     use_r_header = FALSE
                     ){
    file <- gsub("\\\\","/",file);
    if (missing(write_header) || write_header == "IfEmpty"){
        write_header <- -1;
    }
    if (ident_from_filename){
        ## Put ident from the filename
        grammar_ident = (gsub("[.][^.]$","",basename(file)));
    }
    if (missing(outputFile)){
        outputFile <- file.path(dirname(file),paste0(basename(file),".d_parser.",write_extension));
    } else if (dir.exists(outputFile)){
        outputFile <- file.path(outputFile,paste0(basename(file),".d_parser.",write_extension));
    }
    outputFile <- gsub("\\\\","/",outputFile);
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
          as.integer(token_type),
          as.integer(use_r_header));
    return(invisible());
}

updateDparser <- function(){ # nocov start
    if (!file.exists(devtools::package_file("src/dparser"))){
        owd <- getwd();
        setwd(devtools::package_file("src"));
        system("git clone --depth 1 https://github.com/jplevyak/dparser")
        setwd(owd);
    }
    if (file.exists(devtools::package_file("src/dparser"))){

        owd <- getwd();
        setwd(devtools::package_file("src/dparser"));
        build <- gsub("([^ ]*) .*","\\1",system("git show-ref heads/master",intern=TRUE));
        setwd(owd);
        cat(sprintf("dparser build %s\n", build));
        for (f in c("dparse.h", "dparse_tables.h", "dsymtab.h", "gram.h", "gramgram.h", "lex.h",
                    "lr.h", "mkdparse.h", "parse.h", "read_binary.h", "scan.h", "util.h", "write_tables.h",
                    "gram.c", "grammar.g.c", "lex.c", "lex.h", "lr.c", "mkdparse.c", "scan.c", "symtab.c",
                    "util.c", "write_tables.c", "d.h")){
            cat(sprintf("\tf: %s\n", f));
            unlink(devtools::package_file("src/", f));
            d <- readLines(devtools::package_file("src/dparser/", f));
            if (f == "d.h"){
                w <- which(regexpr('#define (REALLOC|MALLOC|FREE|CALLOC)', d) != -1);
                d <- d[-w];
                w <- which(regexpr('#include "arg.h"', d) != -1);
                d <- d[-w];
                w <- which(regexpr('#include <string.h>', d) != -1)
                d <- c(d[1:w], "#include <R.h>", "#include <Rinternals.h>",
                       "#define D_MAJOR_VERSION 1",
                       "#define D_MINOR_VERSION 30",
                       sprintf('#define D_BUILD_VERSION "R-%s"', build),
                       d[seq(1 + w, length(d))]);
            }
            d <- gsub("REALLOC", "R_chk_realloc", d);
            d <- gsub("MALLOC[(]", "R_chk_calloc(1,", d);
            d <- gsub("CALLOC[(]", "R_chk_calloc(1,", d);
            d <- gsub("FREE", "Free", d);
            d <- gsub("([ \t])printf[(]", "\\1Rprintf(", d);
            if (f == "write_tables.c"){
                w  <- which(regexpr('#include "dparse_tables.h"', d, fixed=TRUE) != -1);
                d <- c(d[1:w], "int d_use_r_headers = 0;", d[seq(w + 1, length(d))]);
                w <- which(regexpr("Available at http://dparser.sf.net", d, fixed=TRUE) != -1);
                for (i in rev(w)){
                    fp <- gsub(" *fprintf[(] *([^ ]*) *,.*", "\\1", d[i]);
                    d[i] = sprintf('fprintf(%s, "\t/*\\n  Generated by RxODE\\\'s mkdparse a port of Make DParser Version %%s\\n", ver);
\tfprintf(%s,"  RxODE available at https://github.com/hallowkm/RxODE\\n");
\tfprintf(%s,"  Original dparser Available at http://dparser.sf.net\\n*/\\n\\n\\n");
\tif (d_use_r_headers){
\t\tfprintf(%s,"\\n\\n#include <R.h>\\n#include <Rinternals.h>\\n#define printf Rprintf\\n\\n");
\t}',fp, fp, fp, fp);
                    d <- d[-(i - 1)];
                }
            }
            if (f == "util.c"){
                w <- which(regexpr("^ *d_warn *[(]", d) != -1);
                w2 <- w + 1;
                while (regexpr("}", d[w2]) == -1){
                    w2 <- w2 + 1;
                }
                d  <- c(d[seq(1, w - 1)],
                        'd_warn(const char *str, ...) {
  char nstr[256];
  char outstr[256*2];
  va_list ap;
  va_start(ap, str);
  snprintf(nstr, 255, "%s", str);
  vsprintf(outstr, nstr, ap);
  va_end(ap);
  warning(outstr);
}',
d[seq(w2 + 1, length(d))]);
                w <- which(regexpr("^ *d_fail *[(] *const +char", d) != -1);
                w2 <- w + 1;
                while (regexpr("}", d[w2]) == -1){
                    w2 <- w2 + 1;
                }
                d  <- c(d[seq(1, w - 1)],
                        'd_fail(const char *str, ...) {
  char nstr[256];
  char outstr[256*2];
  va_list ap;
  va_start(ap, str);
  snprintf(nstr, 255, "Parser Fail: %s", str);
  vsprintf(outstr, nstr, ap);
  va_end(ap);
  error(outstr);
}',
d[seq(w2 + 1, length(d))]);
            }

            sink(devtools::package_file("src/", f));
            cat(paste(d, collapse="\n"));
            sink()
        }
    }
} # nocov end

updatePkg <- function(){ # nocov start
    cat("Update README\n");
    owd <- getwd();
    on.exit({setwd(owd)});
    setwd(devtools::package_file());
    knitr::knit(devtools::package_file("README.Rmd"))
    cat("Update Parser c file\n");
    mkdparse(devtools::package_file("inst/tran.g"),
             devtools::package_file("src/"),
             grammar_ident="RxODE");
    file <- gsub("^([#]line [0-9]+ )\".*(src)/+(.*)\"","\\1\"\\2/\\3\"",
                 readLines(devtools::package_file("src/tran.g.d_parser.c")))
    sink(devtools::package_file("src/tran.g.d_parser.c"))
    cat(paste(file,collapse="\n"));
    sink();
    unlink(devtools::package_file("src/tran.o"))
    sink(devtools::package_file("R/version.R"))
    cat("##\' Version and repository for this RxODE package.
##\'
##\' @return A character vector with the version and repository.
##\' @author Matthew L. Fidler
##\' @keywords internal
##\' @export
rxVersion <- function(){return(c(version=\"");
    cat(gsub("Version: +", "", readLines(devtools::package_file("DESCRIPTION"), 2)[2]))
    cat("\",repo=\"");
    cat("https://github.com/")
    tmp <- readLines(devtools::package_file(".git/config"))
    cat(gsub("\\.git$", "", gsub(".*git@github.com:", "", tmp[which(tmp == '[remote "origin"]')[1]+1])))
    cat("\"))}\n");
    sink();
    devtools::load_all();
} # nocov end
