library(shiny)
library(shinythemes)

# Load the 5-fold xgboost models and recalibration models for severe case among UK Biobank population
load(paste0("cohortC_reCaliModels.Rdata"))
severe_cali.models <- cali.models
severe_xgb.models <- xgb.models

# Load the 5-fold xgboost models and recalibration models for fatal case among UK Biobank population
load(paste0("cohortD_reCaliModels.Rdata"))
fatal_cali.models <- cali.models
fatal_xgb.models <- xgb.models

# Load all the predicted and calibrated risk among UK Biobank population
severe.sample.risk <- read.csv("cohortC_clinical_risk.beta_recal.tsv", header = TRUE, 
                               stringsAsFactors=FALSE, sep="\t")

fatal.sample.risk <- read.csv("cohortD_clinical_risk.beta_recal.tsv", header = TRUE, 
                              stringsAsFactors=FALSE, sep="\t")

# Load xgboost explainer model (Fold 1 only)
load("cohortC_explainer.RData")
severe.xgb <- xgb
severe.explainer <- explainer

load("cohortD_explainer.RData")
fatal.xgb <- xgb
fatal.explainer <- explainer

# Load the feature name lookup table
df.nlookup <- read.csv("names_lookup.csv", header = TRUE, stringsAsFactors=FALSE, quote="\"")


# Define UI for dataset viewer app ----
ui <- fluidPage(theme = shinytheme("journal"), # For detail, pls refer to https://bootswatch.com/journal/

  # App title ----
  titlePanel("Risk Estimation for Covid-19"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Include clarifying text ----
      h4(class="text-info", "Please input your clinical risk factors",
               "and press the button to estimate the Relative Risk",
               "and Ranking of predicted risk if infected by Covid-19."),
      
      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change). This is useful if
      # the computations required to render output are inordinately
      # time-consuming.
      actionButton("estimate", "Estimate Outcome", class = "btn btn-info"),

      hr(),
      
      # Input: Specify Sex ----
      radioButtons("sex", "Gender:",
                  choices = c("male"=1, "female"=0),
                  inline=TRUE),

      # Input: Specify current Age ----
      #tags$div(title="Range: 49 - 86",
      #  numericInput("age", "Age (Years), valid range is Age 49 to 86:", 0, min=49, max=86)
      #),
      sliderInput("age", "Age (Years):",
                  min = 49, max = 86,
                  value = 0),

      tags$div(title="Please click here to slide through",
               sliderInput("waist_circum", "Waist circumference (cm):",
                  min = 20, max = 197,
                  value = 90)
      ),
      sliderInput("weight", "Weight (kg):",
                  min = 5, max = 200,
                  value = 78),
      sliderInput("height", "Height (cm):",
                  min = 50, max = 200,
                  value = 170),
      sliderInput("dia_blood_pressure", "Diastolic Blood pressure:",
                  min = 32, max = 148,
                  value = 82),
      sliderInput("sys_blood_pressure", "systolic Blood pressure:",
                  min = 62, max = 268,
                  value = 138),
      
      radioButtons("ethnic", "Ethnic background:",
                   choices = c("White" = 1,
                               "Mixed" = 2,
                               "Asian/Asian British" = 3,
                               "Black/Black British" = 4,
                               "Chinese" = 5,
                               "Others" = 6),
                   selected=1,
                   inline=TRUE),      
      
      radioButtons("alcohol_status", "Alcohol drinker status:",
                   choices = c("Never" = 0,
                               "Previous" = "1",
                               "Current" = "2"),
                   selected=0,
                   inline=TRUE),
      
      radioButtons("smoking_status", "Smoking status:",
                   choices = c("Never" = 0,
                               "Previous" = "1",
                               "Current" = "2"),
                   selected=0,
                   inline=TRUE),
      
      h4("Diagnosis by Doctor:"),
      checkboxInput("dx_cancer", "Cancer", FALSE),
      checkboxInput("AF", "Atrial fibrillation (AF)", FALSE),
      checkboxInput("CAD", "Coronary artery disease (CAD)", FALSE),
      checkboxInput("Cognitive_impairment", "Cognitive_impairment", FALSE),
      checkboxInput("Depression", "Depression", FALSE),
      checkboxInput("Haemorr_stroke", "Hemorrhagic stroke", FALSE),
      checkboxInput("Anxiety", "Anxiety disorder", FALSE),
      checkboxInput("Psychosis", "Psychosis", FALSE),
      checkboxInput("Bipolar", "Bipolar disorder", FALSE),
      checkboxInput("Heart_failure", "Heart Failure", FALSE),
      checkboxInput("Asthma_Doc_Self", "Asthma", FALSE),
      checkboxInput("COPD_Doc_Self", "Chronic obstructive pulmonary disease (COPD)", FALSE),
      checkboxInput("T1DM_ICD_Self", "Type 1 diabetes mellitus (T1DM)", FALSE),
      checkboxInput("T2DM_ICD_Self", "Type 2 diabetes mellitus (T2DM)", FALSE),

      sliderInput("cnt_cancers", "Number of cancers:", min=0, max=6, value=0),
      sliderInput("cnt_noncancer", "Number of non-cancer illnesses:", min=0, max=29, value=0),
      sliderInput("cnt_treatment", "Number of tratments/medications taken:", min=0, max=36, value=0),
     
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(type="tabs",
                  tabPanel("Result",
                           list(
                             # Output: Table of relative risk and ranking result
                             withTags(
                               div(class="jumbotron",
                                   h4(class="lead", "Risk estimation in the population"),
                                   h4(tableOutput("result"))
                                  )
                               ),

                             # Output: The predicted relative risk and the ranking ----

                             withTags(
                               div(class="jumbotron",
                                     h4(class="lead", "The following refers to the predicted risks in the population"),
                                     h4("Risk of being infected and developing a severe disease (hospitalized/fatal):"),
                                     h4(style="padding-left:2em", htmlOutput("summingup1")),
                                     h4("Risk of being infected and developing a fatal disease:"),
                                     h4(style="padding-left:2em", htmlOutput("summingup2"))
                                   )
                               ),
                             
                           
                             # Paper reference
                             withTags(
                               div(class="jumbotron", checked=NA,
                                   h4(class="lead", 'Please refer to following for more details.'),
                                   h4(a("Uncovering clinical risk factors and prediction of severe COVID-19: A machine learning approach based on UK Biobank data",
                                                          href = "https://doi.org/10.1101/2020.09.18.20197319"),
                                     "2020 Kenneth C.Y. Wong, Yong Xiang, Hon-Cheong So 10.1101/2020.09.18.20197319")
                                   )
                               ),
                           
                             # Medical Disclaimer
                             withTags(
                               div(class="jumbotron", checked=NA,
                                   h4(class="lead", strong("Medicial Disclaimer:")),
                                   h4("This tool is designed for general educational and research purposes only and is not intended in any way ",
                                     "to substitute for professional medical advice, consultation, diagnosis, or treatment. Any analysis, report, ",
                                     "or information contained in or produced by this tool is intended to serve as a supplement to, and not a ",
                                     "substitute for the knowledge, expertise, skill and judgment of health care professionals. The contents of ",
                                     "this tool may be of interest to medical professions or other health care providers. Medical professionals ",
                                     "and health care providers should exercise their own judgment in determining whether this tool is appropriate ",
                                     'for their practice or their patients. The information accessed through this tool is provided "AS IS" and ',
                                     "without any warranties, whether expressed or implied. In no event shall this tool under this Agreement, be ",
                                     "considered to be in any form, medical care, treatment, or therapy for patients or users of this tool. ",
                                     "The creators of this tool shall not be liable for any direct, indirect, consequential, special, exemplary, ", 
                                     "lost profits, or other damages incurred by the user of this tool.")
                                   )
                               )
                           )
                  ),
                  tabPanel("Shapley Plot",
                           list(
                             withTags(
                               div(class="jumbotron",
                                     h4(class="lead",  "Contribution of clinical factors to the risk of being infected and developing a severe disease."),
                                     plotOutput("severe_shap_dc_plot")
                                   )
                               ),
                             withTags(
                               div(class="jumbotron",
                                     h4(class="lead", "Contribution of clinical factors to the risk of being infected and developing a fatal disease."),
                                     plotOutput("fatal_shap_dc_plot")
                                   )
                               )
                             )
                           )

      )
    )
  )
)




# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  library(xgboost)
  library(betacal)
  library(xgboostExplainer)
  library(ggplot2)
  
  modelC_path <- "cohortC_reCaliModels.Rdata"
  modelD_path <- "cohortD_reCaliModels.Rdata"
  
  # Prevalence in the UK biobank sample (updated at 14 Dec 2020)
  severe_preval <- 0.005097
  fatal_preval <- 0.0010242
  
  # Return the requested dataset ----
  # Note that we use eventReactive() here, which depends on
  # input$update (the action button), so that the output is only
  # updated when the user clicks the button
  getEstimation <- eventReactive(input$estimate, {
    print(paste0("Call getEstimate(): ", input$age))
    if (input$age > 0) {
      # Prepare the test dataset
      df.test <- data.frame(
        sex.1=input$sex,
        alcohol_status.1=ifelse(input$alcohol_status==1, 1, 0),
        alcohol_status.2=ifelse(input$alcohol_status==2, 1, 0),
        dx_cancer.1=ifelse(input$dx_cancer==TRUE, 1, 0),
        ethnic.2=ifelse(input$ethnic==2, 1, 0),
        ethnic.3=ifelse(input$ethnic==3, 1, 0),
        ethnic.4=ifelse(input$ethnic==4, 1, 0),
        ethnic.5=ifelse(input$ethnic==5, 1, 0),
        ethnic.6=ifelse(input$ethnic==6, 1, 0),
        smoking_status.1=ifelse(input$smoking_status==1, 1, 0),
        smoking_status.2=ifelse(input$smoking_status==2, 1, 0),
        AF.1=ifelse(input$AF==TRUE, 1, 0),
        CAD.1=ifelse(input$CAD==TRUE, 1, 0),
        Cognitive_impairment.1=ifelse(input$Cognitive_impairment==TRUE, 1, 0),
        Depression.1=ifelse(input$Depression==TRUE, 1, 0),
        Haemorr_stroke.1=ifelse(input$Haemorr_stroke==TRUE, 1, 0),
        Anxiety.1=ifelse(input$Anxiety==TRUE, 1, 0),
        Psychosis.1=ifelse(input$Psychosis==TRUE, 1, 0),
        Bipolar.1=ifelse(input$Bipolar==TRUE, 1, 0),
        Heart_failure.1=ifelse(input$Heart_failure==TRUE, 1, 0),
        Asthma_Doc_Self.1=ifelse(input$Asthma_Doc_Self==TRUE, 1, 0),
        COPD_Doc_Self.1=ifelse(input$COPD_Doc_Self==TRUE, 1, 0),
        T1DM_ICD_Self.1=ifelse(input$T1DM_ICD_Self==TRUE, 1, 0),
        T2DM_ICD_Self.1=ifelse(input$T2DM_ICD_Self==TRUE, 1, 0),
        current_age=input$age,
        cnt_cancers=input$cnt_cancers,
        cnt_noncancer=input$cnt_noncancer,
        cnt_treatment=input$cnt_treatment,
        waist_circum=input$waist_circum,
        weight=input$weight,
        bmi=input$weight/((input$height/100.0)^2),
        dia_blood_pressure=input$dia_blood_pressure,
        sys_blood_pressure=input$sys_blood_pressure
      )
      print(df.test)
      df.test.num <- data.matrix(df.test)

      # Re-calibrate and measure on test set
      N <- 5
      severe_recal_probs <- c()
      fatal_recal_probs <- c()
      for (i in 1:N) {
        # Pre-calibrated predicted probability of test dataset
        severe_test_res <- predict(severe_xgb.models[[i]], df.test.num)
        fatal_test_res <- predict(fatal_xgb.models[[i]], df.test.num)
        
        # Calibrated predicted probability of test dataset
        severe_recal_prob = beta_predict(severe_test_res, severe_cali.models[[i]])
        fatal_recal_prob = beta_predict(fatal_test_res, fatal_cali.models[[i]])
        
        # Save the predicted and re-calibrated risk 
        severe_recal_probs <- c(severe_recal_probs, severe_recal_prob)
        print(paste0("Severe Risk recalibration: ", severe_test_res, " -> ", severe_recal_prob))
        fatal_recal_probs <- c(fatal_recal_probs, fatal_recal_prob)
        print(paste0("Fatal Risk recalibration: ", fatal_test_res, " -> ", fatal_recal_prob))
      }
      
      avg.severe_recal_probs <- mean(severe_recal_probs)
      severe_rr <- avg.severe_recal_probs/severe_preval
      avg.fatal_recal_probs <- mean(fatal_recal_probs)
      fatal_rr <- avg.fatal_recal_probs/fatal_preval

      print(paste0("Avg severe recal. prob: ", round(avg.severe_recal_probs,5)))
      print(paste0("Avg fatal recal. prob: ", round(avg.fatal_recal_probs,5)))
      
      
      # Evaluate the risk ranking of the test subject
      severe_test_data.risk <- data.frame(eid="9999999", clinical_risk=avg.severe_recal_probs, outcome="ctrl")
      fatal_test_data.risk <- data.frame(eid="9999999", clinical_risk=avg.fatal_recal_probs, outcome="ctrl")
      
      df.risk <- rbind(severe.sample.risk, severe_test_data.risk)
      df.risk$rank <- rank(df.risk$clinical_risk)
      position <- df.risk[df.risk$eid==9999999,4]
      severe_top_rank_in_percent <- 100*(1- position/nrow(df.risk))
      
      df.risk <- rbind(fatal.sample.risk, fatal_test_data.risk)
      df.risk$rank <- rank(df.risk$clinical_risk)
      position <- df.risk[df.risk$eid==9999999,4]
      fatal_top_rank_in_percent <- 100*(1- position/nrow(df.risk))

    } else {
      severe_rr <- 0
      severe_top_rank_in_percent <- 0
      fatal_rr <- 0
      fatal_top_rank_in_percent <- 0
    }
    rr <- c(severe_rr, fatal_rr)
    rank <- c(severe_top_rank_in_percent, fatal_top_rank_in_percent)
    df.result <- data.frame(rr, rank)
    colnames(df.result) <- c("Relative Risk", "Rank in top k %")
    rownames(df.result) <- c("severe case", "fatal case")
    df.result
  }, ignoreNULL = FALSE)
  
  
  getTestData <- eventReactive(input$estimate, {
    print(paste0("Call getTestData(): ", input$age))
    if (input$age > 0) {
      # Prepare the test dataset
      df.test <- data.frame(
        sex.1=input$sex,
        alcohol_status.1=ifelse(input$alcohol_status==1, 1, 0),
        alcohol_status.2=ifelse(input$alcohol_status==2, 1, 0),
        dx_cancer.1=ifelse(input$dx_cancer==TRUE, 1, 0),
        ethnic.2=ifelse(input$ethnic==2, 1, 0),
        ethnic.3=ifelse(input$ethnic==3, 1, 0),
        ethnic.4=ifelse(input$ethnic==4, 1, 0),
        ethnic.5=ifelse(input$ethnic==5, 1, 0),
        ethnic.6=ifelse(input$ethnic==6, 1, 0),
        smoking_status.1=ifelse(input$smoking_status==1, 1, 0),
        smoking_status.2=ifelse(input$smoking_status==2, 1, 0),
        AF.1=ifelse(input$AF==TRUE, 1, 0),
        CAD.1=ifelse(input$CAD==TRUE, 1, 0),
        Cognitive_impairment.1=ifelse(input$Cognitive_impairment==TRUE, 1, 0),
        Depression.1=ifelse(input$Depression==TRUE, 1, 0),
        Haemorr_stroke.1=ifelse(input$Haemorr_stroke==TRUE, 1, 0),
        Anxiety.1=ifelse(input$Anxiety==TRUE, 1, 0),
        Psychosis.1=ifelse(input$Psychosis==TRUE, 1, 0),
        Bipolar.1=ifelse(input$Bipolar==TRUE, 1, 0),
        Heart_failure.1=ifelse(input$Heart_failure==TRUE, 1, 0),
        Asthma_Doc_Self.1=ifelse(input$Asthma_Doc_Self==TRUE, 1, 0),
        COPD_Doc_Self.1=ifelse(input$COPD_Doc_Self==TRUE, 1, 0),
        T1DM_ICD_Self.1=ifelse(input$T1DM_ICD_Self==TRUE, 1, 0),
        T2DM_ICD_Self.1=ifelse(input$T2DM_ICD_Self==TRUE, 1, 0),
        current_age=input$age,
        cnt_cancers=input$cnt_cancers,
        cnt_noncancer=input$cnt_noncancer,
        cnt_treatment=input$cnt_treatment,
        waist_circum=input$waist_circum,
        weight=input$weight,
        bmi=input$weight/((input$height/100.0)^2),
        dia_blood_pressure=input$dia_blood_pressure,
        sys_blood_pressure=input$sys_blood_pressure
      )
      data.matrix(df.test)
    } else {
      NULL
    }
    
  }, ignoreNULL = FALSE)
  
  output$result <- renderTable({
    getEstimation()
  }, rownames = TRUE)
  
  output$summingup1 <- renderText({
    df <- getEstimation()
    print(df)
    severe_rr <- sprintf("%.1f", df[1,1])
    severe_top_rank_in_percent <- sprintf("%.1f", df[1,2])
    severe_top_rank_in_1kppl <- sprintf("%.0f", df[1,2]*10)
    paste0('<b>You are ranked at the <font color=\"#FF0000\">', severe_top_rank_in_percent, "</font>",
           "th percentile (i.e. ranked ", '<font color=\"#FF0000\">', severe_top_rank_in_1kppl, 
           '</font> among <font color=\"#FF0000\">1000</font> people, ',
           'ordering from the highest to lowest risk) and have <font color=\"#FF0000\">', severe_rr, 
           "</font> times the mean risk among UK Biobank (UKBB) participants.</b>")
  })
  
  output$summingup2 <- renderText({
    df <- getEstimation()
    fatal_rr <- sprintf("%.1f", df[2,1])
    fatal_top_rank_in_percent <- sprintf("%.1f", df[2,2])
    fatal_top_rank_in_1kppl <- sprintf("%.0f", df[2,2]*10)
    paste0('<b>You are ranked at the <font color=\"#FF0000\">', fatal_top_rank_in_percent, "</font>",
           "th percentile (i.e. ranked ", '<font color=\"#FF0000\">', fatal_top_rank_in_1kppl, 
           '</font> among <font color=\"#FF0000\">1000</font> people, ',
           'ordering from the highest to lowest risk) and have <font color=\"#FF0000\">', fatal_rr, 
           "</font> times the mean risk among UK Biobank (UKBB) participants.</b>")
  })

  output$severe_shap_dc_plot <- renderPlot({
    df.test.num <- getTestData()
    if (!is.null(df.test.num)) {
      
      # Prepare the arguments for waterfall plot
      dtest.rename <- xgb.DMatrix(df.test.num)
      df.test.num.rename <- df.test.num
      
      # Rename the predictor names
      colnames(dtest.rename) <- df.nlookup$desc[match(colnames(dtest.rename), df.nlookup$predictor)]
      colnames(df.test.num.rename) <- df.nlookup$desc[match(colnames(df.test.num), df.nlookup$predictor)]
      
      pred.breakdown = explainPredictions(severe.xgb, severe.explainer, df.test.num.rename)
      
      df.test.num.rename <- apply(df.test.num.rename, c(1, 2), round, 4)
      showWaterfall(severe.xgb, severe.explainer, dtest.rename, data.matrix(df.test.num.rename) , 1, type="binary", threshold=0.03)+
        theme(axis.text.x = element_text(angle = -70, hjust = 0))
    } else {
      NULL
    }

  })
  
  output$fatal_shap_dc_plot <- renderPlot({
    df.test.num <- getTestData()
    if (!is.null(df.test.num)) {
      
      # Prepare the arguments for waterfall plot
      dtest.rename <- xgb.DMatrix(df.test.num)
      df.test.num.rename <- df.test.num
      
      # Rename the predictor names
      colnames(dtest.rename) <- df.nlookup$desc[match(colnames(dtest.rename), df.nlookup$predictor)]
      colnames(df.test.num.rename) <- df.nlookup$desc[match(colnames(df.test.num), df.nlookup$predictor)]
      
      pred.breakdown = explainPredictions(fatal.xgb, fatal.explainer, df.test.num.rename)
      
      df.test.num.rename <- apply(df.test.num.rename, c(1, 2), round, 4)
      showWaterfall(fatal.xgb, fatal.explainer, dtest.rename, data.matrix(df.test.num.rename) , 1, type="binary", threshold=0.03)+
        theme(axis.text.x = element_text(angle = -70, hjust = 0))
    } else {
      NULL
    }
    
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)

## How to deploy an application to http://labso.shinyapps.io
# 1. install.packages('rsconnect')
# 2. rsconnect::setAccountInfo(name='labso', token='81AEBC3E32E67956FF5A49F49B1DEF60', secret='<SECRET>')
# 3. library(rsconnect)
# 4. rsconnect::deployApp('/home/kenneth/RA/covid19/Rshiny/covid19_risk_pred')
# 5. Visit https://labso.shinyapps.io/covid19_risk_pred/