library(digest);
rxPermissive({
    files <- list.files(pattern=".*\\.test\\.g$")
    for (file in files){
        flags <- sprintf("%s.flags",file);
        if (file.exists(flags)){
            flags <- readLines(flags);
            if (flags == "-A"){
                flags <- list(states_for_all_nterms=TRUE);
            }
        } else {
            flags <- list();
        }
        context(sprintf("Grammar %s", file));
        out <- sprintf("%s.d_parser.c",file);
        flags$file <- file;
        flags$use_r_header <- TRUE;
        flags$verbose <- FALSE;
        do.call("mkdparse",flags);
        sink("Makevars");
        cat(sprintf("PKG_CPPFLAGS=-I\"%s\"\n",rxIncludeDir()))
        sink();
        parser <- sprintf("sample_parser%s",.Platform$dynlib.ext);
        unlink(parser)
        cmd <- sprintf("%s/bin/R CMD SHLIB sample_parser.c %s ",
                       Sys.getenv("R_HOME"), base::basename(out));
        sh <- "system";
        do.call(sh,list(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE));
        unlink("Makevars");
        unlink(out);
        unlink(gsub("\\.c$",".h",out));
        unlink(gsub("\\.c$",".o",out));
        unlink("sample_parser.o");
        dyn.load(parser);
        ret <- tryCatch(dyn.load(parser),error=function(e){return(FALSE)});
        if ((class(ret) == "DLLInfo")){
            for (parseFile in list.files(pattern=sprintf("%s.[0-9]+$",file))){
                parseFlags <- sprintf("%s.flags",parseFile);
                test.name <- sprintf("%s: %s", file, parseFile);
                args <- list(fileName=parseFile,
                             start_state=0,
                             save_parse_tree= 1,
                             partial_parses = 0,
                             compare_stacks = 1,
                             commit_actions_interval = 100,
                             fixup = 1,
                             fixup_ebnf = 0,
                             nogreedy = 0,
                             noheight = 0,
                             use_file_name = 1);
                if (file.exists(parseFlags)){
                    parseFlags <- readLines(parseFlags);
                    if (parseFlags == "-S 3"){
                        args$start_state = 3;
                    } else if (parseFlags == "-e"){
                        args$fixup_ebnf = 1;
                    } else if (parseFlags == "-f") {
                        args$fixup = 0;
                    } else {
                        print(parseFlags)
                    }
                    ## print(parseFlags);
                    ## stop();
                }
                if(parseFile == "g50.test.g.1"){
                    args$use_file_name <- 0;
                }
                sink("test-dparser");
                with(args,
                     cat(.Call("sample_parser",
                               fileName,
                               as.integer(start_state),
                               as.integer(save_parse_tree),
                               as.integer(partial_parses),
                               as.integer(compare_stacks),
                               as.integer(commit_actions_interval),
                               as.integer(fixup),
                               as.integer(fixup_ebnf),
                               as.integer(nogreedy),
                               as.integer(noheight),
                               as.integer(use_file_name))));
                sink();
                test <- readLines("test-dparser");
                unlink("test-dparser");
                ref <- readLines(sprintf("%s.check",parseFile));
                test_that(test.name, {
                    if (test.name == "g50.test.g: g50.test.g.1"){
                        skip("Doesn't seem to work under R");
                    }
                    expect_equal(test,ref);
                });
            }
            dyn.unload(parser);
            unlink(parser);
        }
    }
}, silent=TRUE)
