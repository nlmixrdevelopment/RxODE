##' Run Shiny app on RxODE object
##'
##' @param object  RxODE family of objects.
##' @param ... Other arguments passed to rxShiny
##' @return Nothing; Starts a shiny server
##' @author Zufar Mulyukov and Matthew L. Fidler
rxShiny <- function(object, ...){
    UseMethod("rxShiny");
}

##' @rdname rxShiny
rxShiny.rxSolve <- function(object, ...){
    ui<- eval(bquote(shinyUI(fluidPage(
        tags$style(HTML("input:invalid {background-color: #FFCCCC;}")),
        tags$script('
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
  Shiny.onInputChange("EnterPressed", Math.random());
  }});
  '
),
fluidRow(
    column(12,
           fluidRow(
               column(5,
                      textAreaInput(
                          inputId = "ode",
                          label = "ODE",
                          width = 360,
                          height = 200,
                          resize = 'both',
                          value = .(rxNorm(object))),

actionButton("goButton", "Compile", align='right'),
hr(),
h4('Time sampling'),
tags$div(style="display: inline-block;",textInput('tmin', 'start', width=60, value = 0)),
tags$div(style="display: inline-block;",textInput('tmax', 'end'  , width=60, value = 100)),
tags$div(style="display: inline-block;",textInput('step', 'step' , width=60, value = 1)),
hr(),
h4('Dosing'),
tags$div(style="display: inline-block;",textInput('dose', 'amount', width=60, value=1)),
tags$div(style="display: inline-block;",textInput('rate', 'rate'  , width=60)),
div(style="display: inline-block;",uiOutput("dosing_cmt")),
br(),
tags$div(style="display: inline-block;",textInput('start'   , 'start'    , width=60, value = 0)),
tags$div(style="display: inline-block;",textInput('ndoses' ,  'repeat'   , width=60, value = 1)),
tags$div(style="display: inline-block;",textInput('interval', 'interval' , width=60, value = 0))

),
column(width = 7,
       fluidRow(
           h4('Initial values'),
           uiOutput("cmts"),
           br(),
           h4('Parameters'),
           uiOutput("pars"),
           uiOutput('message', placeholder = FALSE),
           ## So may parens;  looks like lisp :)
           uiOutput("plotTabs"))))))))))

    server <- function(input, output, session) {
        values <- reactiveValues()
        values$m1 = NULL
        values$res = NULL
        values$msg = 'NULL'


        tmp=tempfile()

        observeEvent(input$goButton, {
            values$m1=NULL
            values$res <- NULL

            values$msg <- capture.output(
                tryCatch({

                    values$m1<-RxODE(model = input$ode, wd = tmp)
                    values$cmts = values$m1$get.modelVars()$state
                    values$pars = values$m1$get.modelVars()$params

                }, error = function(e) return(e))
            )

        })
        observe({if(!is.null(values$m1)) values$msg = 'Model is compiled. Enter parameters.'})

        output$message <- renderPrint({

            if(length(values$msg)==0)
                return(invisible())
            if( values$msg=='NULL' )
                return(invisible())

            cat('<h3>Message:</h3>')
            cat(paste(gsub('<|>','',values$msg), collapse = '<br>'))
        })

        output$dosing_cmt <- renderUI({
            if (is.null(values$cmts))
                return()
            selectInput("into", "into", width=80, choices = values$cmts, selected=values$cmts[1])
        })

        output$cmts <- renderUI({
            if (is.null(values$cmts))
                return()
            if(length(values$cmts)==0)
                return('none')

            lapply(values$cmts, function(x) {
                div(style="display: inline-block;", textInput(x, x, value = 0, width=60))
            })
        })

        output$pars <- renderUI({
            if (is.null(values$pars))
                return()
            if(length(values$pars)==0)
                return('none')

            lapply(values$pars, function(x) {
                div(style="display: inline-block;", textInput(x, x, value = 1, width=60))
            })
        })

        output$plotTabs = renderUI({
            if (is.null(values$res))
                return(tabsetPanel())
            dat=values$res

            cmts =colnames(values$res)[-1]

            sel.tab = isolate(ifelse(is.null(input$plot.tabs),cmts[1], input$plot.tabs))

            tabs =
                lapply(cmts, function(cmt){
                    plotname <- paste("plot", cmt, sep="")
                    output[[plotname]] <- renderPlot({
                        ggplot(as.data.frame(dat), aes_(x=as.name("time"), y=as.name(cmt))) + geom_line()
                    })

                    tabPanel(cmt,  plotOutput(plotname))
                })
            do.call(tabsetPanel, c(tabs, id='plot.tabs', selected = sel.tab))
        })



        observeEvent({input$EnterPressed},{
            values$msg <- capture.output(
                tryCatch({
                    solveODE()
                }, error = function(e) return(e))
            )
        })


        solveODE <- function(){
            values$res <- NULL
            if (is.null(values$m1))
                return()

            stime = as.numeric(input$tmin)
            etime = as.numeric(input$tmax)
            tstep = as.numeric(input$step)
            dose  = as.numeric(input$dose)
            rate  = as.numeric(input$rate)

            if(is.na(rate)){rate <- NULL}
            if(is.na(dose)){dose=0}

            ndoses   = as.numeric(input$ndoses)
            start    = as.numeric(input$start)
            interval = as.numeric(input$interval)

            into = ifelse(!is.null(input$into),1,match(input$into, values$cmts))

            ev <- eventTable()

            ev$add.sampling(seq(stime, etime, tstep ))

            ev$add.dosing(dose = dose, start.time=start,
                          nbr.doses= ndoses, rate=rate,
                          dosing.interval = interval,
                          dosing.to = into)

            params <- .(object$params.single)

            cmts = values$cmts

            inits <- NULL
            if(length(cmts)>0){
                init_str=paste0(cmts, '=as.numeric(input$',cmts,')')
                init_str = paste('c(',toString(init_str),')')
                inits <- eval(parse(text=init_str))
            }

            pars = values$pars
            params = NULL
            if(length(pars)>0){
                param_str = paste0(pars,'=as.numeric(input$',pars,')')
                param_str = paste('c(',toString(param_str),')')
                params <- eval(parse(text=param_str))
            }

            values$res =
                values$m1$solve(params, ev, inits)
        }

        session$onSessionEnded(stopApp)
    }


    shinyApp(ui = ui, server = server)
}


