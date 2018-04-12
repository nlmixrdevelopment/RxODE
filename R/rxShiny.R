##' Use Shiny to help develop an RxODE model
##'
##' @param object A RxODE family of objects. If not supplied a
##'     2-compartment indirect effect model is used.  If it is
##'     supplied, use the model associated with the RxODE object for
##'     the model exploration.
##' @param ... Other arguments passed to rxShiny.  Currently doesn't
##'     do anything.
##' @return Nothing; Starts a shiny server
##' @author Zufar Mulyukov and Matthew L. Fidler
rxShiny <- function(object, ...){
    UseMethod("rxShiny");
}
##' @rdname rxShiny
rxShiny.default <- function(object="
C2 = centr/V2
C3 = peri/V3
d/dt(depot) =-KA*depot
d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3
d/dt(peri)  =                    Q*C2 - Q*C3
d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff
", ...){
    ui<- eval(bquote(shiny::shinyUI(shiny::fluidPage(
                                               shiny::tags$style(shiny::HTML("input:invalid {background-color: #FFCCCC;}")),
        shiny::tags$script('
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
  shiny::Shiny.onInputChange("EnterPressed", Math.random());
  }});
  '
),
shiny::fluidRow(
    shiny::column(12,
           shiny::fluidRow(
               shiny::column(5,
                      shiny::textAreaInput(
                          inputId = "ode",
                          label = "ODE",
                          width = 360,
                          height = 200,
                          resize = 'both',
                          value=.(rxNorm(object))),

shiny::actionButton("goButton", "Compile", align='right'),
shiny::hr(),
shiny::h4('Time sampling'),
shiny::tags$div(style="display: inline-block;",shiny::textInput('tmin', 'start', width=60, value = 0)),
shiny::tags$div(style="display: inline-block;",shiny::textInput('tmax', 'end'  , width=60, value = 100)),
shiny::tags$div(style="display: inline-block;",shiny::textInput('step', 'step' , width=60, value = 1)),
shiny::hr(),
shiny::h4('Dosing'),
shiny::tags$div(style="display: inline-block;",shiny::textInput('dose', 'amount', width=60, value=1)),
shiny::tags$div(style="display: inline-block;",shiny::textInput('rate', 'rate'  , width=60)),
shiny::div(style="display: inline-block;",shiny::uiOutput("dosing_cmt")),
shiny::br(),
shiny::tags$div(style="display: inline-block;",shiny::textInput('start'   , 'start'    , width=60, value = 0)),
shiny::tags$div(style="display: inline-block;",shiny::textInput('ndoses' ,  'repeat'   , width=60, value = 1)),
shiny::tags$div(style="display: inline-block;",shiny::textInput('interval', 'interval' , width=60, value = 0))

),
shiny::column(width = 7,
       shiny::fluidRow(
           shiny::h4('Initial values'),
           shiny::uiOutput("cmts"),
           shiny::br(),
           shiny::h4('Parameters'),
           shiny::uiOutput("pars"),
           shiny::uiOutput('message', placeholder = FALSE),
           shiny::uiOutput("plotTabs")
       )))))))))

    server <- function(input, output, session) {
        values <- reactiveValues()
        values$m1 = NULL
        values$res = NULL
        values$msg = 'NULL'


        tmp=tempfile()

        shiny::observeEvent(input$goButton, {
            values$m1=NULL
            values$res <- NULL

            values$msg <- capture.output(
                tryCatch({

                    values$m1<-RxODE(model = input$ode, wd = tmp)
                    values$cmts = values$m1$get.modelVars()$state
                    values$pars = rxParams(values$m1)

                }, error = function(e) return(e))
            )


        })
        observe({if(!is.null(values$m1)) values$msg = 'Model is compiled. Enter parameters.'})

        output$message <- shiny::renderPrint({

            if(length(values$msg)==0)
                return(invisible())
            if( values$msg=='NULL' )
                return(invisible())

            cat('<h3>Message:</h3>')
            cat(paste(gsub('<|>','',values$msg), collapse = '<br>'))
        })

        output$dosing_cmt <- shiny::renderUI({
            if (is.null(values$cmts))
                return()
            selectInput("into", "into", width=80, choices = values$cmts, selected=values$cmts[1])
        })

        output$cmts <- shiny::renderUI({
            if (is.null(values$cmts))
                return()
            if(length(values$cmts)==0)
                return('none')

            lapply(values$cmts, function(x) {
                shiny::div(style="display: inline-block;", shiny::textInput(x, x, value = 0, width=60))
            })
        })

        output$pars <- shiny::renderUI({
            if (is.null(values$pars))
                return()
            if(length(values$pars)==0)
                return('none')

            lapply(values$pars, function(x) {
                shiny::div(style="display: inline-block;", shiny::textInput(x, x, value = 1, width=60))
            })
        })

        output$plotTabs = shiny::renderUI({
            if (is.null(values$res))
                return(tabsetPanel())
            dat=values$res

            cmts =colnames(values$res)[-1]

            sel.tab = isolate(ifelse(is.null(input$plot.tabs),cmts[1], input$plot.tabs))

            tabs =
                lapply(cmts, function(cmt){
                    plotname <- paste("plot", cmt, sep="")
                    output[[plotname]] <- renderPlot({
                        ggplot2::ggplot(as.data.frame(dat), ggplot2::aes_(x=as.name("time"), y=as.name(cmt))) + ggplot2::geom_line()
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

            params <-
                c(KA = .291, CL = 18.6,
                  V2 = 40.2, Q = 10.5, V3 = 297.0,
                  Kin = 1.0, Kout = 1.0, EC50 = 200.0)

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


    shiny::shinyApp(ui = ui, server = server)
}
