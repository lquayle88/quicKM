## quicKM web app
## Lewis Quayle
## 21/01/2022


### GLOBAL

## load libraries

library("bslib")
library("cowplot")
library("dplyr")
library("ggplot2")
library("readr")
library("shiny")
library("shinyWidgets")
library("stringr")
library("survival")
library("survminer")
library("summaryBox")
library("tidyr")


## load data

load(file.path("www", "tcga_brca.RData"))


## define utility functions

# filter counts data

filter_count_data <- function(clinical.data, expr.data, annotation, genes.input) {
  
  # extract vector of case ids from filtered clinical data
  
  case.ids <- pull(clinical.data, id)
  
  # switch annotation type
  
  annotation <- switch(annotation,
                       "Ensembl ID" = "ensembl",
                       "Entrez ID" = "entrez",
                       "HUGO Symbol" = "symbol")
  
  # create a character vector used to filter count data rows
  
  if (genes.input == "") {
    
    genes <- NULL
    
  } else {
    
    genes <- genes.input %>%
      strsplit(split = paste(",", ";", "\\s", sep = "|")) %>%
      unlist() %>%
      toupper()
    
  }
  
  # filter gene level expression data
  
  expr.data %>%
    dplyr::select(ensembl, entrez, symbol, all_of(case.ids)) %>%
    dplyr::filter(get(annotation) %in% genes)
  
}

# calculate geometric mean

calc_geo_mean <- function(x) {
  
  x %>%
    log() %>%
    mean() %>%
    exp()
  
}

# generate survival data frame

generate_survival_data <- function(clinical.data, expr.data, split.cohort, censor.threshold, threshold) {
  
  # only calculate with data frame resulting from valid genes input
  
  req(nrow(expr.data) > 0)
  
  # create an expression matrix
  
  expr.data <- expr.data %>%
    dplyr::select(-c(ensembl, entrez, symbol))
  
  # extract vector of expression values for a single gene or compute the geometric mean for multiple genes
  
  expr.value <- sapply(expr.data, FUN = calc_geo_mean)
  
  # create the data frame used in survival analysis
  
  surv.data <- data.frame("case.id" = names(expr.value),
                          "os_status" = clinical.data$os_status,
                          "os_time" = clinical.data$os_time,
                          "expr.value" = expr.value,
                          row.names = NULL)
  
  # order by increasing gene / signature expression value
  
  surv.data <- arrange(surv.data, expr.value)
  
  # calculate quantiles in expression data
  
  quant <- quantile(x = surv.data$expr.value, probs = seq(from = 0, to = 1, by = 0.25))
  
  # define the quantile cutoff according to user input
  
  cutoff <- switch(split.cohort, "lower quartile" = quant[2], "median" = quant[3], "upper quartile" = quant[4])
  
  # binary stratum classifier variable
  
  surv.data <- surv.data %>%
    mutate(stratum = ifelse(surv.data$expr.value < cutoff, 0, 1))
  
  # define follow-up period
  
  if (censor.threshold) {
    
    surv.data %>%
      mutate(os_status = ifelse(os_time <= threshold, os_status, 0))
    
  } else {
    
    surv.data %>%
      dplyr::filter(os_time <= threshold)
    
  }
  
}

# clinical data plots

plot_clinical_svar <- function(data, x.var, x.lab) {
  
  ggplot(data = data, aes_string(x = x.var)) +
    geom_bar(colour = "#4D4D4D", fill = "#1E90FF") +
    theme_bw() +
    theme(aspect.ratio = 3/4,
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(10, 0, 10, 0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0, 15, 0, 15)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    labs(x = x.lab, y = "Number of Cases")
  
}

plot_clinical_mvar <- function(data, x.var, grp.var, x.lab, grp.lab) {
  
  ggplot(data = data, aes_string(x = x.var, alpha = grp.var)) +
    geom_bar(position = "dodge", colour = "#4D4D4D", fill = "#1E90FF") +
    scale_alpha_manual(values = c(1, 0.6, 0.3)) +
    theme_bw() +
    theme(aspect.ratio = 3/4,
          plot.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 14, face = "bold", margin = margin(10, 0, 10, 0)),
          axis.title.y = element_text(size = 14, face = "bold", margin = margin(0, 15, 0, 15)),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          legend.background = element_rect(colour = "#4D4D4D")) +
    labs(x = x.lab, y = "Number of Cases", alpha = grp.lab)
  
}

# survival analysis plot

plot_survival <- function(data, title, pval, conf.int, censor, surv.median.line, threshold) {
  
  # custom plot theme
  
  theme_survival <- function() {
    
    theme_bw() +
      theme(plot.background = element_blank(),
            plot.title = element_text(size = 18, face = "bold"),
            panel.grid.major = element_line(size = rel(0.5)), 
            panel.grid.minor = element_line(size = rel(0.25)),
            axis.title.x = element_text(size = 14, face = "bold", margin = margin(10, 0, 10, 0)),
            axis.title.y = element_text(size = 14, face = "bold", margin = margin(0, 15, 0, 15)),
            axis.text = element_text(size = 12),
            axis.ticks = element_line(size = rel(0.5)),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.background = element_rect(colour = "#4D4D4D"))
    
  }
  
  # custom risk table theme
  
  theme_rt <- function() {
    
    theme_survival() +
      theme(plot.title = element_text(size = 16, face = "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())
    
  }
  
  # compute survival curves using Kaplan-Meier model formula
  
  fit <- survfit(data = data, formula = Surv(time = os_time, event = os_status) ~ stratum)
  
  # generate plot
  
  srvplot <- ggsurvplot(fit = fit,
                        data = data,
                        pval = pval,
                        pval.size = 4.5,
                        pval.coord = c(1, 0.1),
                        conf.int = conf.int,
                        censor = censor,
                        risk.table = FALSE,
                        surv.median.line = switch(surv.median.line, "TRUE" = "hv", "FALSE" = "none"),
                        ggtheme = theme_survival(),
                        palette = c("#053061", "#E41A1C"),
                        axes.offset = TRUE,
                        break.time.by = 1,
                        xlim = c(0, threshold),
                        ylim = c(0, 1),
                        break.y.by = 0.10,
                        title = title,
                        xlab = "Time (Years)",
                        ylab = "Overall Survival Probability",
                        legend = "right",
                        legend.title = "Group",
                        legend.labs = c("Low", "High"))
  
  srvplot <- srvplot[[1]]
  
  # generate risk table
  
  risktab <- ggsurvplot(fit = fit,
                        data = data,
                        risk.table = TRUE,
                        ggtheme = theme_rt(),
                        fontsize = 4,
                        palette = c("#053061", "#E41A1C"),
                        break.time.by = 1,
                        xlim = c(0, threshold),
                        xlab = "Time (Years)",
                        risk.table.title = "Number at Risk",
                        legend.title = "",
                        legend.labs = c("Low", "High"))
  
  risktab <- risktab[[2]]
  
  # align on grid
  
  plot_grid(srvplot,
            risktab,
            ncol = 1,
            align = "v",
            axis = "lr",
            rel_heights = c(8,2))
  
}


### USER-INTERFACE

ui <- navbarPage(
  
  # app title
  
  titlePanel(title = div(img(src = "logo.png", width = 180)), windowTitle = "quicKM: Speedy Survival Analysis"),
  
  # theme
  
  theme = bs_theme(bootswatch = "cerulean"),
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## HOME PAGE
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  tabPanel(title = "Home",
           icon = icon("home"),
           includeMarkdown(file.path("www", "homepage.md")),
           setBackgroundImage("site_bg.png")),
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## BREAST CANCER
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  navbarMenu(title = "Breast", icon = icon("fas fa-project-diagram"),
             
             
             ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             ## TCGA-BRCA
             ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             tabPanel(title = "TCGA-BRCA", icon = icon("far fa-folder"),
                      
                      tabsetPanel(
                        
                        
                        ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                        ## TCGA-BRCA - CLNICAL DATA
                        ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                        
                        tabPanel(title = "Clinical Data",
                                 
                                 sidebarLayout(
                                   
                                   ## SIDEBAR
                                   
                                   sidebarPanel(
                                     
                                     pickerInput(inputId = "race_tcga_brca",
                                                 label = "Race or Ethnicity:",
                                                 choices = levels(tcga_brca_clinical$race),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "menopause_status_tcga_brca",
                                                 label = "Menopausal Status:",
                                                 choices = levels(tcga_brca_clinical$menopause_status),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "histological_type_tcga_brca",
                                                 label = "Histological Type:",
                                                 choices = levels(tcga_brca_clinical$histological_type),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "pathological_stage_tcga_brca",
                                                 label = "Pathological Stage:",
                                                 choices = levels(tcga_brca_clinical$stage),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "er_status_tcga_brca",
                                                 label = "Oestrogen Receptor (ER) Status:",
                                                 choices = levels(tcga_brca_clinical$er_status),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "pr_status_tcga_brca",
                                                 label = "Progesterone Receptor (PR) Status:",
                                                 choices = levels(tcga_brca_clinical$pr_status),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     pickerInput(inputId = "her2_status_tcga_brca",
                                                 label = "HER2 Status:",
                                                 choices = levels(tcga_brca_clinical$her2_status),
                                                 multiple = TRUE,
                                                 options = list(title = "select one or more")),
                                     
                                     actionButton(inputId = "ab_1_tcga_brca",
                                                  label = "Confirm Selection",
                                                  icon = icon("far fa-check-circle"),
                                                  style = "width:100%"),
                                     
                                     uiOutput("dl_1_tcga_brca")
                                     
                                   ), # end sidebarPanel
                                   
                                   ## MAIN PANEL
                                   
                                   mainPanel(
                                     
                                     uiOutput("brca_summarybox", style = "margin-top:2em"),
                                     
                                     fluidRow(
                                       column(width = 12, style = "margin-top:2em", plotOutput("clinical_plots_tcga_brca"))
                                     )
                                     
                                   ) # end mainPanel
                                   
                                 ) # end sidebarLayout
                        ), # end tabPanel
                        
                        
                        ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                        ## TCGA-BRCA - SURVIVAL ANALYSIS
                        ## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                        
                        tabPanel(title = "Survival Analysis",
                                 
                                 sidebarLayout(
                                   
                                   # SIDEBAR
                                   
                                   sidebarPanel(
                                     
                                     
                                     prettyRadioButtons(
                                       inputId = "annotation_tcga_brca",
                                       label = "Gene Annotation Type:", 
                                       choices = c("Ensembl ID", "Entrez ID", "HUGO Symbol"),
                                       selected = "HUGO Symbol",
                                       inline = TRUE, 
                                       status = "primary",
                                       fill = TRUE
                                     ),
                                     
                                     fluidRow(
                                       column(width = 10,
                                              style = "display:inline-block;vertical-align:top",
                                              textAreaInput(inputId = "gene_filter_tcga_brca",
                                                            label = "Gene or Gene List:",
                                                            rows = 3)),
                                       column(width = 2,
                                              style = "display:inline-block;vertical-align:top",
                                              circleButton(inputId = "genes_help_tcga_brca",
                                                           size = "xs",
                                                           color = "primary",
                                                           icon = icon("question")))
                                     ),
                                     
                                     pickerInput(inputId = "define_cutoff_tcga_brca",
                                                 label = "Split Patient Cohort By:",
                                                 choices = c("lower quartile", "median", "upper quartile"),
                                                 selected = "median"),
                                     
                                     fluidRow(
                                       column(width = 10,
                                              style = "display:inline-block;vertical-align:top",
                                              sliderInput(inputId = "define_time_tcga_brca",
                                                 label = "Follow-Up Threshold:",
                                                 min = 2,
                                                 max = floor(max(tcga_brca_clinical$os_time)),
                                                 step = 1,
                                                 value = floor(max(tcga_brca_clinical$os_time)))),
                                       
                                       column(width = 2,
                                              style = "display:inline-block;vertical-align:top",
                                              circleButton(inputId = "followup_help_tcga_brca",
                                                           size = "xs",
                                                           color = "primary",
                                                           icon = icon("question")))
                                       ),

                                     materialSwitch(inputId = "censor_threshold_tcga_brca",
                                                    label = "Censor at Threshold:", 
                                                    value = TRUE,
                                                    status = "primary"),
                                     
                                     textInput(inputId = "plot_title_tcga_brca",
                                               label = "Plot Title:",
                                               width = "100%"),
                                     
                                     materialSwitch(inputId = "show_conf_tcga_brca",
                                                    label = "Confidence Intervals:", 
                                                    value = FALSE,
                                                    status = "primary"),
                                     
                                     materialSwitch(inputId = "med_surv_line_tcga_brca",
                                                    label = "Median Survival Line:", 
                                                    value = FALSE,
                                                    status = "primary"),
                                     
                                     materialSwitch(inputId = "censor_tcga_brca",
                                                    label = "Censor Points:", 
                                                    value = TRUE,
                                                    status = "primary"),
                                     
                                     materialSwitch(inputId = "show_pval_tcga_brca",
                                                    label = "P-Value:", 
                                                    value = TRUE,
                                                    status = "primary"),
                                     
                                     actionButton(inputId = "ab_2_tcga_brca",
                                                  label = "Draw KM Plot",
                                                  icon = icon("fas fa-chart-line"),
                                                  style = "width:100%"),
                                     
                                     fluidRow(
                                       column(6, uiOutput("dl_2_tcga_brca")),
                                       column(6, uiOutput("dl_3_tcga_brca"))
                                     ),
                                     
                                     uiOutput("dl_4_tcga_brca")
                                     
                                   ), # end sidebarPanel
                                   
                                   # MAIN PANEL
                                   
                                   mainPanel(
                                     
                                     plotOutput("tcga_brca_survival"),
                                     style = "margin-top:1em"
                                     
                                   ) # end mainPanel
                                   
                                 ) # end sidebarLayout
                        ) # end tabPanel
                      ) # end tabsetPanel
             ), # end tabPanel
             
             
             ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             ## STUDY 2
             ## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(                      
                          
                          
                          
                        )
                      )
             ) # end of tabset
             
  ), # end navbarMenu
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## PROSTATE CANCER
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  navbarMenu(title = "Prostate", icon = icon("fas fa-project-diagram"),
             
             # STUDY 1
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(                          
                          
                          #imageOutput("under_development")
                          
                        )
                      )
             ),
             
             # STUDY 2
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(                      
                          
                          
                          
                        )
                      )
             )
  ), # end of prostate studies
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## LUNG CANCER
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  navbarMenu(title = "Lung", icon = icon("fas fa-project-diagram"),
             
             # STUDY 1
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(
                          
                          
                          
                        )
                      )
             ),
             
             # STUDY 2
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(                      
                          
                          
                          
                        )
                      )
             )
  ), # end of lung studies
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## COLON CANCER
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  navbarMenu(title = "Colon", icon = icon("fas fa-project-diagram"),
             
             # STUDY 1
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(
                          
                          
                          
                        )
                      )
             ),
             
             # STUDY 2
             
             tabPanel(title = "UNDER DEVELOPMENT", icon = icon("far fa-folder"),
                      sidebarLayout(
                        sidebarPanel(
                          
                          
                          
                        ),
                        
                        mainPanel(                      
                          
                          
                          
                        )
                      )
             )
  ), # end of colon studies
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## FAQ PAGE
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  tabPanel(title = "FAQ", icon = icon("far fa-question-circle"),
           
           includeMarkdown(file.path("www", "faq.md"))
           
  ),
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## ABOUT PAGE
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  tabPanel(title = "About", icon = icon("fas fa-book"),
           
           includeMarkdown(file.path("www", "about.md"))
           
  )
  
  
) # end navbarPage


### SERVER

server <- function(input, output) {
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## START-UP MODAL
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # showModal(modalDialog(title = "Welcome to quicKM",
  #                       size = "xl",
  #                       fade = TRUE,
  #                       "PLACE HOLDER FOR USER DEMO",
  #                       footer = modalButton("Get Started"),
  #                       easyClose = TRUE))
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## TCGA-BRCA - CLINICAL DATA
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ## filter clinical data
  
  clinical_selection <- eventReactive(input$ab_1_tcga_brca, {
    
    # race
    
    fil_race <- if (is.null(input$race_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(race) %>%
        unique()
    } else {
      input$race_tcga_brca
    }
    
    # menopause status
    
    fil_menopause <- if (is.null(input$menopause_status_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(menopause_status) %>%
        unique()
    } else {
      input$menopause_status_tcga_brca
    }
    
    # histological type
    
    fil_histological <- if (is.null(input$histological_type_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(histological_type) %>%
        unique()
    } else {
      input$histological_type_tcga_brca
    }
    
    # pathological stage
    
    fil_pathological <- if (is.null(input$pathological_stage_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(stage) %>%
        unique()
    } else {
      input$pathological_stage_tcga_brca
    }
    
    # er status
    
    fil_er <- if (is.null(input$er_status_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(er_status) %>%
        unique()
    } else {
      input$er_status_tcga_brca
    }
    
    # pr status
    
    fil_pr <- if (is.null(input$pr_status_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(pr_status) %>%
        unique()
    } else {
      input$pr_status_tcga_brca
    }
    
    # her2 status
    
    fil_her2 <- if (is.null(input$her2_status_tcga_brca)) {
      tcga_brca_clinical %>%
        pull(her2_status) %>%
        unique()
    } else {
      input$her2_status_tcga_brca
    }
    
    # filter clinical data
    
    tcga_brca_clinical %>%
      dplyr::filter(race %in% fil_race,
                    menopause_status %in% fil_menopause,
                    histological_type %in% fil_histological,
                    stage %in% fil_pathological,
                    er_status %in% fil_er,
                    pr_status %in% fil_pr,
                    her2_status %in% fil_her2)
    
  })
  
  
  ## update action button 1 once run
  
  observeEvent(input$ab_1_tcga_brca, {
    
    updateActionButton(inputId = "ab_1_tcga_brca",
                       label = "Update Selection")
    
  })
  
  
  ## clinical data case number summary box
  
  output$brca_summarybox <- renderUI({
    
    data <- clinical_selection()
    case.number <- nrow(data)
    
    fluidRow(
      
      summaryBox(title = "Number of Cases Selected Using Current Parameters:",
                 value = case.number,
                 width = 12,
                 icon = "fas fa-file-medical",
                 style = "primary",
                 border = "left")
      
    )
    
  })
  
  
  ## clinical data plots
  
  # reactive graphic object
  
  clinical_plots_tcga_brca <- reactive({
    
    p1 <- plot_clinical_svar(data = clinical_selection(), x.var = "race", x.lab = "Race")
    p2 <- plot_clinical_svar(data = clinical_selection(), x.var = "menopause_status", x.lab = "Menopause Status")
    p3 <- plot_clinical_svar(data = clinical_selection(), x.var = "histological_type", x.lab = "Histological Type")
    p4 <- plot_clinical_svar(data = clinical_selection(), x.var = "stage", x.lab = "Pathological Stage")
    
    input.data <- clinical_selection() %>%
      dplyr::select(id, er_status, pr_status, her2_status) %>%
      pivot_longer(cols = ends_with("status"),names_to = "receptor", values_to = "status") %>%
      mutate(receptor = toupper(str_remove(receptor, "_status")),
             receptor = factor(receptor, levels = c("ER", "PR", "HER2"), ordered = TRUE))
    
    
    p5 <- plot_clinical_mvar(data = input.data,
                             x.var = "receptor",
                             grp.var = "status",
                             x.lab = "Receptor Type",
                             grp.lab = "Expression Status")
    
    plot_grid(p1, p2, p3, p4, p5,
              ncol = 2,
              nrow = 3,
              align = "v",
              axis = "lr")
    
  })
  
  # render graphic object
  
  output$clinical_plots_tcga_brca <- renderPlot({
    
    print(clinical_plots_tcga_brca())
    
  },
  
  bg = "transparent",
  width = 920,
  height = 920
  
  )
  
  
  ## clinical data plots download
  
  output$download_clinical_plots_tcga_brca <- downloadHandler(
    
    filename = function() {
      paste0("tcga_brca_clinical_summary_plot.png")
    },
    
    content = function(file) {
      
      ggsave(file,
             plot = clinical_plots_tcga_brca(),
             scale = 1.5,
             width = 8.27,
             height = 11.69,
             units = "in",
             dpi = 600,
             bg = "#FFFFFF",
             device = "png")
      
    }
    
  )
  
  output$dl_1_tcga_brca <- renderUI({
    
    req(clinical_plots_tcga_brca())
    downloadButton(outputId = "download_clinical_plots_tcga_brca",
                   label = "Download Plots",
                   icon = icon("download"),
                   style = "width:100%;margin-top:1em")
    
  })
  
  
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ## TCGA-BRCA - SURVIVAL ANALYSIS
  ## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ## input gene(s) of interest and filter gene level expression data
  
  counts_selection <- eventReactive(input$ab_2_tcga_brca, {
    
    filter_count_data(clinical.data = clinical_selection(),
                      expr.data = tcga_brca_expression,
                      annotation = input$annotation_tcga_brca,
                      genes.input = input$gene_filter_tcga_brca)
    
  })
  
  
  ## update action button 2 once run
  
  observeEvent(input$ab_2_tcga_brca, {
    
    if (nrow(counts_selection()) == 0) {
      
      showModal(modalDialog(title = 'Input Error: "Gene or Gene List"',
                            includeMarkdown(file.path("www", "error_gene_input.md")),
                            easyClose = TRUE))
      
    } else {
      
      updateActionButton(inputId = "ab_2_tcga_brca",
                         label = "Update Parameters")
      
    }
    
  })
  
  
  ## gene selection help modal
  
  observeEvent(input$genes_help_tcga_brca, {
    
    showModal(modalDialog(title = "quicKM Help",
                          includeMarkdown(file.path("www", "help_gene_input.md")),
                          easyClose = TRUE)
    )
    
  })

  ## follow-up threshold help modal
  
  observeEvent(input$followup_help_tcga_brca, {
    
    showModal(modalDialog(title = "quicKM Help",
                          includeMarkdown(file.path("www", "help_followup.md")),
                          easyClose = TRUE)
    )
    
  })
  
  
  ## create survival analysis data frame
  
  surv_data <- eventReactive(input$ab_2_tcga_brca, {
    
    generate_survival_data(clinical.data = clinical_selection(),
                           expr.data = counts_selection(),
                           split.cohort = isolate(input$define_cutoff_tcga_brca),
                           censor.threshold = isolate(input$censor_threshold_tcga_brca),
                           threshold = isolate(input$define_time_tcga_brca))
    
  })
  
  
  ## survival analysis plot
  
  # reactive graphic object
  
  surv_plot_tcga_brca <- reactive({
    
    plot_survival(data = surv_data(),
                  title = isolate(input$plot_title_tcga_brca), 
                  pval = isolate(input$show_pval_tcga_brca),
                  conf.int = isolate(input$show_conf_tcga_brca),
                  censor = isolate(input$censor_tcga_brca),
                  surv.median.line = isolate(input$med_surv_line_tcga_brca),
                  threshold = isolate(input$define_time_tcga_brca))
    
  })
  
  # render graphic object
  
  output$tcga_brca_survival <- renderPlot({
    
    print(surv_plot_tcga_brca())
    
  },
  
  bg = "transparent",
  height = 700
  
  )
  
  
  ## expression data download
  
  output$download_expr_data_tcga_brca <- downloadHandler(
    
    filename = function() {
      paste0("tcga_brca_input_gene_expression_data.csv")
    },
    
    content = function(file) {
      readr::write_csv(x = counts_selection(), file)
    }
    
  )
  
  output$dl_2_tcga_brca <- renderUI({
    
    req(surv_data())
    downloadButton(outputId = "download_expr_data_tcga_brca",
                   label = "Expression Data",
                   icon = icon("download"),
                   style = "width:100%;margin-top:1em")
    
  })
  
  
  ## survival analysis data download
  
  output$download_surv_data_tcga_brca <- downloadHandler(
    
    filename = function() {
      paste0("tcga_brca_survival_analysis_data.csv")
    },
    
    content = function(file) {
      write_csv(x = surv_data(), file)
    }
    
  )
  
  output$dl_3_tcga_brca <- renderUI({
    
    req(surv_data())
    downloadButton(outputId = "download_surv_data_tcga_brca",
                   label = "Survival Data",
                   icon = icon("download"),
                   style = "width:100%;margin-top:1em")
    
  })
  
  
  ## survival analysis KM plot download
  
  output$download_surv_plot_tcga_brca <- downloadHandler(
    
    filename = function() {
      
      paste0("km_plot.png")
      
    },
    
    content = function(file) {
      
      ggsave(file,
             plot = surv_plot_tcga_brca(),
             device = "png",
             width = 11.69,
             height = 8.27,
             units = "in",
             dpi = 600,
             bg = "#FFFFFF")
      
    }
    
  )
  
  output$dl_4_tcga_brca <- renderUI({
    
    req(surv_plot_tcga_brca())
    downloadButton(outputId = "download_surv_plot_tcga_brca",
                   label = "KM Plot",
                   icon = icon("download"),
                   style = "width:100%;margin-top:1em")
    
  })
  
  
} # end of server


## RUN APP

shinyApp(ui, server)

