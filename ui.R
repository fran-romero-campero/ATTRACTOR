# Define UI for application that draws a histogram
ui <- fluidPage(
  ##shinythemes::themeSelector(),
  theme = shinytheme("sandstone"),
  ##theme = shinytheme("simplex"),
  ##theme = shinytheme("journal"),
  
  # Application title, introductory text and sidebar navigation
  fluidRow(
    column(
      width = 2,
      img(src='attractor_logo_2.jpg', align = "center", width=200),
      tags$br(),
      radioButtons(inputId = "navigation_bar", width="100%",selected="home",
                   label="",
                   choices=c(
                     "Home" = "home",
                     "Individual gene analysis" = "individual_gene",
                     "Multiple transcription factor analysis" = "multiple_gene",
                     "Tutorials" = "tutorials",
                     "GitHub repository" = "github",
                     "Citation" = "citation"
                   ))),
    column(
      width = 8,
      tags$div(align = "center", 
               tags$h1(tags$b("ATTRACTOR,"),tags$i("Arabidopsis Thaliana"), "TRanscriptionAl Circadian neTwORk")),
      tags$br(),tags$br(),
      
      conditionalPanel(condition = "input.navigation_bar == 'home'",
                       tags$div(align="justify", "The", tags$b("circadian clock"), "and", tags$b("light signalling"), "play central roles in", 
                                tags$i("plant physiology"), "and", tags$i("development."), "As a consequence, massive amounts of",
                                tags$b("omics data"), " have been generated to characterize their individual components. Nonetheless, 
                                these data remain fragmented and researchers who want to explore the synergistic regulation exerted by the 
                                circadian clock and light signalling need to consult different papers and resources making imperative 
                                the use of", tags$b("molecular systems biology"), "techniques to integrate and make easily accesible 
                                all the generated information."),
                       tags$div(align="justify", tags$b("ATTRACTOR,"),"is a web based tool for the analysis of the synergistic transcriptional control 
                                exerted by the circadian clock and light signalling over genes exhibiting rhythmic expression profiles in the model plant ", 
                                tags$i(tags$b("Arabidopsis thaliana.")), tags$b("ATTRACTOR,"), "consists of a ", tags$b("transcriptional network"), 
                                " that integrates transcriptomic data collected over diurnal cycles with 12 hours of light and 12 hours of darkness 
                                with cistromic data generated using ChIP-seq for key transcriptional factors and regulators in the circadian clock 
                                and light signalling. Specifically, our network is composed of 5778 nodes or genes with diurnal rhythmic expression profiles and
                                14529 edges or transcriptional regulations. The transcription factors and regulators included in our network comprise the
                                components of the morning and central loops", tags$b("CCA1, LHY,"), "the pseudo response regulator family members", 
                                tags$b("TOC1,PRR5, PRR7"), "and ", tags$b("PRR9;"), "as well as some components of the evening loop such as", tags$b("LUX, ELF3"), "and", tags$b("ELF4."),
                                "In order to capture synergistic regulations with light signalling we added the light sensors and transcriptional regulators phytochromes",
                                tags$b("PHYA"), "and", tags$b("PHYB,"),"the cryptochrome ", tags$b("CRY2"), "as well as the light transcriptional factors from the phytochrome interacting factor
          family", tags$b("PIF5, PIF4"), " and ", tags$b("PIF3."),"Finally, the phytochrome interacting transcriptional factor", tags$b("FHY1"), "(Far-red elongated Hypocotyl 1)
          is also included in our network."),
                       tags$div(align="justify","Use the navigation bar on the left to explore the different utilities in ATTRACTOR or alternatively",
                                tags$a(href="https://www.youtube.com/watch?v=8eJN5zrMZbI", target="_blank", tags$b("view our video tutorial."))),
                       tags$br(), tags$br()#,
                       # actionButton("run", "Run Animation"),
                       # fluidRow(column(width = 9,
                       #                 plotOutput("networkAnimation")),
                       #          column(width = 3,
                       #                 tags$br(), tags$br(),tags$br(), tags$br(),tags$br(), tags$br(),tags$br(), tags$br(),
                       #                 plotOutput("clockAnimation")))
                       # plotOutput("networkAnimation"),
                       # tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(), tags$br(),
                       # plotOutput("clockAnimation")
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'github'",
                       tags$div(align = "justify", tags$b("ATTRACTOR,"), "is entirely developed using 
        the R package", tags$b( tags$a(href="https://shiny.rstudio.com/", "shiny.")), "The 
        source code is released under", tags$b("GNU General Public License v3.0"), "and is hosted at",
                                tags$b("GitHub."), "If you experience any problem using ATTRACTOR please create an", tags$b(tags$a(href="https://github.com/fran-romero-campero/ATTRACTOR/issues","issue")), "in GitHub and we will address it."),
                       tags$div(align="center",tags$h1(tags$b(tags$a(href="https://github.com/fran-romero-campero/ATTRACTOR", "ATTRACTOR at GitHub"))))
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'citation'",
                       tags$div(align = "justify", "We are strongly committed to", tags$b("open access software"), 
                                "and", tags$b("open science."),"Following our philosophy we have deposited our GitHub code 
                       into", tags$a(href="https://zenodo.org/record/3780022#.XqwzNvlS9uQ", target="_blank",tags$b("Zenodo")), ", a
                       general-purpose open-access repository developed under the", 
                                tags$a(href="https://www.openaire.eu/", target="_blank", tags$b("European OpenAIRE program.")), "Meanwhile we publish 
                       our work in a journal if you find", tags$b("ATTRACTOR"), "useful in your research we would be most grateful if you cite 
                       our GitHub repository with a,", tags$b("DOI"),  "as follows:",
                                tags$br(),
                                tags$br(),
                                tags$div(tags$h4(tags$b("de los Reyes, P., Romero-Losada, A.B., Romero-Campero, F.J. (2020) ATTRACTOR, 
                       Arabidopsis Thaliana TRanscriptionAl Circadian neTwORk v1.0, Zenodo, doi:10.5381/zenodo.3780022")))),
                       
                       tags$br(),
                       tags$br(),
                       #tags$div(align="center", img(src='smiley.png', align = "center", width=200,hight=200)),
                       tags$br()
                       
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'tutorials'",
                       tags$div(align="center",uiOutput("video_tutorial")),
                       tags$div(align = "justify", 
                                tags$br(),
                                tags$br(),
                                tags$div(tags$h4(tags$b("Above you can find a video tutorial on how to use the different tools implemented 
                                in ATTRACTOR to explore the transcriptional regulation exerted by the circadian clock and light signalling in 
                                the model plant", tags$i("Arabidopsis thaliana")))))
                       
      ),
      
      conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
                       tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                 transcription factors or regulators over an", tags$b("individually selected"), "gene as well as the effect observed in its
                 expression profile. Follow the steps below:",
                                tags$ol(
                                  tags$li("Select a specific gene from our network using the", tags$b("Target Gene"),
                                          "dropdown menu on the left below. You can enter either the AGI identifier or primary symbol for the",
                                          "gene of interest."), 
                                  tags$li("Select several transcription factors (TFs) to explore their regulation over the
                           previously selected target gene using the",tags$b("Select Transcription Factors"), "checkboxes on the left below."), 
                                  tags$li("Results will be depicted on the tabs below. The", tags$b("Clock Visualizer tab"), "shows a circular 
                            representation of the diurnal cycle with the selected TFs and target gene located at
                            the time point of their expression peak. Binding of TFs to the target gene and 
                            the observed effect over its expression are represented using colored arrows. The", tags$b("Expression 
                            Visualizer tab"),  "shows the target gene expression profile in a linear representation with the selected
                            TFs located over the time point of their expression peak. The binding of 
                            the selected TFs to the target gene promoter and the observed effect is represented using colored arrows.
                            The", tags$b("Peak Visualizer tab"), "shows the selected target gene genomic location and the peaks detected
                            in our analysis of the correspondig TFs ChIP-seq data. Here you can specify the length of the gene promoter and 
                            3' region as well as DNA TFs binding motifs to search for in the detected peak regions with the specified score
                            for identity.")
                                  # tags$li("On the tab", tags$b("Clock Visualizer"), "a circular representation of the diurnal cycle under study
                                  #         will be depicted. The selected transcription factors will be located close to the edge of the circle at the
                                  #         specific time point when the corresponding ChIP-seq data were generated. The selected target gene will be 
                                  #         located towards the center of the circle pointing at the specific time point where its expression is highest. 
                                  #         An arrow will be drawn from a transcription factor to the target gene when the corresponding transcription factor
                                  #         binds to the promoter of the target gene. The node of the transcription factors and arrows will be colored red to represent a decrease in the target gene expression
                                  #         after binding of the transcription factor, green to represent an increase or grey when no change is observed.")
                                ))
      ),
      conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                       tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to explore the coordinated regulation of several 
                                transcription factors or regulators over their common targets. Follow the steps below:", 
                                
                                tags$ol(
                                  tags$li("Select your TFs of interest using the", tags$b("Select Transcription Factors"),
                                          "checkbox menu on the left below."),
                                  tags$li("You can also select a cluster of genes exhibiting a specific rhythmic pattern of expression by chosing
                                          the time point for their peak and trough using the dropdown menus", 
                                          tags$b("Select a specific rhythmic gene expression pattern with peak at: ... and trough at:")),
                                  tags$li("Check the box", tags$b("Visualize Edges"), "when you want to depict arrows from TFs to their target genes."),
                                  tags$li("Click on the", tags$b("SELECT GENES"), "to visualize your selected TFs common target genes exhibiting the 
                                           specified rhythmic expression pattern in our transcriptional network. Explore the different tabs to 
                                           download a table with the selected genes, perform a signficance analysis of the overlap between the selected 
                                           TFs targets and the specified expression pattern as well as", tags$b("GO term, pathways 
                                           and binding motifs enrichment analysis"), ".")
                                )
                       )
      )
      
    ),
    column(
      width = 2,
      img(src='logo_ibvf.jpg', align = "center", width=100),
      img(src='logo_us.png', align = "center", width=100),
      tags$br(),tags$br(),tags$br(),
      img(src='logo_csic.jpg', align = "center", width=100),
      tags$br(),
      tags$br(),
      tags$div(align = "center", width=60,
               #HTML("<script type=\"text/javascript\" id=\"clstr_globe\" src=\"//cdn.clustrmaps.com/globe.js?d=_vn8-gT1zeKYdG0FCV7CsZkw0zisYiGXkw3ZoGonziU\"></script>")
               #HTML("<script type=\"text/javascript\" src=\"//rf.revolvermaps.com/0/0/1.js?i=5lpcyfvr1gd&amp;s=220&amp;m=0&amp;v=false&amp;r=false&amp;b=000000&amp;n=false&amp;c=ff0000\" async=\"async\"></script>")
               HTML("<script type=\"text/javascript\" src=\"//rf.revolvermaps.com/0/0/2.js?i=5eb04j6ugth&amp;m=0&amp;s=130&amp;c=ff0000&amp;t=1\" async=\"async\"></script>")
      )
    )
  ),
  
  ## Separation for different tools
  tags$br(),tags$br(),
  
  ## Conditional panel for individual gene analysis
  conditionalPanel(condition = "input.navigation_bar == 'individual_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Select target gene to study
                            selectizeInput(inputId = "target.gene",
                                           label = "Target Gene:",
                                           choices = genes.selectize,
                                           multiple = FALSE),
                            ## Check box for the TF peaks to represent
                            checkboxGroupInput(inputId = "selected.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = list("PHYA ZT00", "PHYB ZT00", "ELF3 ZT00", "CCA1 ZT02", "LHY ZT02", "ELF3 ZT04","FHY1 ZT04", "PIF4 ZT04", "PIF5 ZT04", "PRR9 ZT04", "CRY2 ZT08", "PIF3 ZT08", "ELF4 ZT10", "PRR5 ZT10", "LUX ZT10", "PRR7 ZT12", "LUX ZT12","CCA1 ZT14", "TOC1 ZT15"),
                                               inline = TRUE,width = "100%")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Clock Visualizer",
                                                 tags$div(align="center",
                                                          plotOutput(outputId = "clock",width = 600, height=600))),
                                        tabPanel(title = "Expression Visualizer", 
                                                 tags$div(align="center",
                                                          plotOutput(outputId = "expression",width = 600, height=600))),
                                        tabPanel(title = "Peak Visualizer",
                                                 column(wellPanel(
                                                   ## Numeric input for promoter length
                                                   numericInput(inputId = "promoter.length",
                                                                label = "Promoter Length",
                                                                value = 2000,
                                                                min = 500,
                                                                max = 2000,
                                                                step = 100),
                                                   ## Numeric input for 3' length
                                                   numericInput(inputId = "threeprime.length",
                                                                label = "3' Length",
                                                                value = 500,
                                                                min = 100,
                                                                max = 500,
                                                                step = 100)),width=3),
                                                 column(wellPanel(
                                                   ## Selectize to choose target gene to represent
                                                   selectizeInput(inputId = "selected.motifs",
                                                                  label = "Select Motifs",
                                                                  selected ="EE",
                                                                  choices = motif.names,
                                                                  multiple = TRUE),
                                                   
                                                   ## Checkbox to select all available motifs
                                                   checkboxInput(inputId = "all.motifs",
                                                                 label = "Select All Motifs:",
                                                                 value = FALSE),
                                                   
                                                   ## Numeric input for PWM score
                                                   numericInput(inputId = "min.score.pwm", 
                                                                label = "Motif Identification Score:",
                                                                value = 100, 
                                                                min = 80,
                                                                max = 100,
                                                                step = 5)),width=9),
                                                 
                                                 #actionButton(inputId = "go",label = "GO"),
                                                 
                                                 fluidRow(
                                                   column(
                                                     plotOutput(outputId = "peak_plot"),
                                                     width=12)
                                                 ))
                            )
                     )
                   )
  ),
  
  
  ## Conditional panel for multiple transcription factors and regulators analysis
  conditionalPanel(condition = "input.navigation_bar == 'multiple_gene'",
                   fluidRow(
                     column(width = 3,
                            ## Check box for the TFs to analyse
                            checkboxGroupInput(inputId = "selected.multiple.tfs",
                                               label = "Select Transcription Factors:",
                                               choices = list("PHYA ZT00", "PHYB ZT00", "ELF3 ZT00", "CCA1 ZT02", 
                                                              "LHY ZT02", "ELF3 ZT04","FHY1 ZT04", "PIF4 ZT04", 
                                                              "PIF5 ZT04", "PRR9 ZT04", "CRY2 ZT08", "PIF3 ZT08", 
                                                              "ELF4 ZT10", "PRR5 ZT10", "LUX ZT10", "PRR7 ZT12", 
                                                              "LUX ZT12","CCA1 ZT14", "TOC1 ZT15"),
                                               inline = TRUE,width = "100%"),
                            tags$b("Select a specific rhythmic gene expression pattern with peak at:"),
                            selectInput(inputId = "peak", label="", 
                                        choices = c("Any ZT" = "any", "ZT0" = "peak0", "ZT4" = "peak4", "ZT8" = "peak8",
                                                    "ZT12" = "peak12", "ZT16" = "peak16", "ZT20" = "peak20"),
                                        selected = NULL, multiple = FALSE, selectize = TRUE),
                            selectInput(inputId = "trough", label="and trough at:", 
                                        choices = c("Any ZT" = "any", "ZT0" = "trough0", "ZT4" = "trough4", "ZT8" = "trough8",
                                                    "ZT12" = "trough12", "ZT16" = "trough16", "ZT20" = "trough20"),
                                        selected = NULL,
                                        multiple = FALSE, selectize = TRUE),
                            checkboxInput(inputId =  "edges",label = "Visualize Edges",value = FALSE),
                            actionButton(inputId = "go_multiple",label = "Select Genes")
                     ),
                     
                     column(width = 9,
                            tabsetPanel(type = "tabs",
                                        tabPanel(title = "Network Visualization",
                                                 plotOutput("networkPlot")
                                        ),
                                        tabPanel(title = "Gene Table",
                                                 tags$br(),
                                                 tags$br(),
                                                 textOutput("empty_overlap_message_1"),
                                                 dataTableOutput(outputId = "outputTable"),
                                                 uiOutput(outputId = "download_ui_for_table")
                                        ),
                                        tabPanel(title = "Overlap Significance",
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$div(align="justify", "In this section, we present the results of a significance analysis of the
                                                           overlap between the targets of the selected transcription factors and gene
                                                           with a specific expresion pattern."),
                                                 tags$br(),
                                                 textOutput("empty_overlap_message_2"),
                                                 textOutput("overlap.message"),
                                                 textOutput("overlap.significance.text"),
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$div(align="center",
                                                          plotOutput("venn.diagram.plot"))
                                        ),
                                        tabPanel(title = "Functional Enrichment",
                                                 tabsetPanel(type = "tabs",
                                                             tabPanel(title = "GO Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a GO term
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      textOutput("empty_overlap_message_3"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "go.background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )), tags$br(),
                                                                      actionButton(inputId = "goterm",label = "GO terms analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      textOutput(outputId = "no_go_results"),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div',h3('Please be patient, computing GO enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "GO table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "textGOTable"),
                                                                                           tags$br(), tags$br(),
                                                                                           dataTableOutput(outputId = "output_go_table"),
                                                                                           htmlOutput(outputId = "revigo"),
                                                                                           uiOutput(outputId = "download_ui_for_go_table")
                                                                                  ),
                                                                                  tabPanel(title = "GO map",
                                                                                           # withSpinner(ui_element =
                                                                                           # plotOutput(outputId = "gomap"),type = 4)),
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "gomap_text"),
                                                                                           tags$br(),
                                                                                           div(style= "overflow:scroll; height:500px;",
                                                                                               addSpinner(plotOutput("gomap"), spin = "circle", color = "#E41A1C"))),
                                                                                  tabPanel(title = "GO barplot",
                                                                                           tags$br(),
                                                                                           htmlOutput(outputId = "barplot_text"),
                                                                                           tags$br(),
                                                                                           # withSpinner(ui_element = 
                                                                                           #   plotOutput(outputId = "bar.plot"),type = 4)),#,inline=TRUE))),
                                                                                           addSpinner(plotOutput("bar.plot"), spin = "circle", color = "#E41A1C")),
                                                                                  tabPanel(title = "GO Enrichment Map",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "emapplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "emap.plot",inline=TRUE))),
                                                                                  tabPanel(title = "GO concept network",
                                                                                           tags$br(), 
                                                                                           htmlOutput(outputId = "cnetplot_text"),
                                                                                           tags$br(),
                                                                                           addSpinner(plotOutput(outputId = "cnet.plot",inline=TRUE)))
                                                                      )
                                                             ),
                                                             tabPanel(title = "KEGG Pathway Enrichment",
                                                                      tags$br(),
                                                                      tags$div(align="justify", "In this section you can perform a KEGG pathways and modules
                                                                               enrichment analysis over the selected genes. First
                                                                               of all, you need to choose the background set of genes between
                                                                               the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                                      tags$br(),
                                                                      textOutput("empty_overlap_message_4"),
                                                                      tags$br(),
                                                                      radioButtons(inputId = "pathway_background", width="100%",selected="allgenome",
                                                                                   label="",
                                                                                   choices=c(
                                                                                     "Complete genome" = "allgenome",
                                                                                     "Genes in network" = "onlynet"
                                                                                   )),
                                                                      actionButton(inputId = "pathway_button",label = "KEGG pathway analysis"),
                                                                      tags$br(),
                                                                      tags$br(),
                                                                      shinyjs::useShinyjs(),
                                                                      hidden(div(id='loading.div.kegg',h3('Please be patient, computing KEGG pathway enrichment ...'))),
                                                                      tags$br(),
                                                                      tabsetPanel(type = "tabs",
                                                                                  tabPanel(title = "Enriched Pathway Table",
                                                                                           tags$br(), tags$br(),
                                                                                           htmlOutput(outputId = "no_kegg_enrichment"),
                                                                                           dataTableOutput(outputId = "output_pathway_table"),
                                                                                           uiOutput(outputId = "download_ui_for_kegg_table")
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Pathway Visualization",
                                                                                           br(),
                                                                                           textOutput("no_pathway_visualization"),
                                                                                           uiOutput(outputId = "kegg_selectize"),
                                                                                           imageOutput("kegg_image"),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
                                                                                           br(), br(), br(), br(), br()
                                                                                  ),
                                                                                  tabPanel(title = "Enriched Module Table",
                                                                                           br(),
                                                                                           htmlOutput(outputId = "text_module_kegg"),
                                                                                           br(), br(),
                                                                                           dataTableOutput(outputId = "output_module_table")
                                                                                  )
                                                                      )
                                                             )
                                                 )
                                        ),
                                        tabPanel(title = "TFBS Enrichment",
                                                 tags$br(),
                                                 tags$div(align="justify", "In this section you can perform a 
                                                          Transcription Factor Binding Sites (TFBS) enrichment analysis 
                                                          over the promoters of the selected genes. First of all, you 
                                                          need to set the some required parameters: the background, 
                                                          the promoter length around the transcriptional start site (TSS) and 
                                                          the motif identification score."),
                                                 tags$br(),
                                                 textOutput("empty_overlap_message_5"),
                                                 tags$br(),
                                                 column(wellPanel(
                                                   ## Select the background
                                                   tags$div(align="justify", "Please, choose the background set of genes between
                                                          the entire genome of", tags$i("Arabidopsis thaliana"), "or just the genes in ATTRACTOR:"),
                                                   # tags$br(),
                                                   radioButtons(inputId = "tfbs_background", width="100%",selected="onlynet",
                                                                label="",
                                                                choices=c(
                                                                  "Complete genome" = "allgenome",
                                                                  "Genes in network" = "onlynet")
                                                   )),width = 4),
                                                 column(wellPanel(
                                                   tags$div(align="justify", "Length of the region upstream of the TSS that will be
                                                              considered gene promoter:"),
                                                   radioButtons(inputId = "up_promoter", width="100%",selected="2000",
                                                                label="",
                                                                choices=c(
                                                                  "500 bp" = "500",
                                                                  "1000 bp" = "1000",
                                                                  "1500 bp" = "1500",
                                                                  "2000 bp" = "2000")
                                                   )),width = 3),
                                                 column(wellPanel(
                                                   tags$div(align="justify", "Length of the region downstream of the TSS that will be
                                                            considered gene promoter:"),
                                                   radioButtons(inputId = "down_promoter", width="100%",selected="500",
                                                                label="",
                                                                choices=c(
                                                                  "0 bp" = "0",
                                                                  "200 bp" = "200",
                                                                  "500 bp" = "500")
                                                   )),width = 3),
                                                 column(wellPanel(
                                                   tags$div(align="justify", "Motif Identification Score:"),
                                                   radioButtons(inputId = "score", width="100%",selected="90",
                                                                label="",
                                                                choices=c(
                                                                  "80" = "80",
                                                                  "85" = "85",
                                                                  "90" = "90",
                                                                  "95" = "95")
                                                   )),width = 2),
                                                 tags$br(),tags$br(),
                                                 tags$br(),tags$br(),tags$br(),
                                                 tags$br(),
                                                 tags$br(),
                                                 tags$br(),
                                                 actionButton(inputId = "tfbs_button",label = "TFBS enrichment analysis"),
                                                 tags$br(),
                                                 textOutput(outputId = "no_tfbs_enrichment"),
                                                 tags$br(),
                                                 shinyjs::useShinyjs(),
                                                 hidden(div(id='loading.div.tfbs',h3('Please be patient, computing TFBS enrichment ...'))),
                                                 tags$br(),
                                                 dataTableOutput(outputId = "output_tfbs_table"),
                                                 uiOutput(outputId = "download_ui_tfbs_table")
                                        ),
                                        tabPanel(title = "Data Retrieval",
                                                 br(), br(),
                                                 tags$div(align="justify", tags$b("ATTRACTOR"), "allows researchers to download bulk data for the 
                                                          set of selected TFs. Users can download for a given list of TFs their gene targets, 
                                                          genomic ChIP signal (BigWig format) or binding sites genomic locations (BED format).
                                                          Mark in the checkbox below the type of data to download, press the RETRIEVE button
                                                          and download the selected data from the generated links:"),
                                                 checkboxGroupInput(inputId = "select_data_retrieval",
                                                                    label = "",
                                                                    choices = list("Target Gene List" = "target_gene_list",
                                                                                   "Genomic Location of Binding Sites (BED format)" = "genomic_locations",
                                                                                   "Genomic Signal (BigWig format)" = "genomic_signal",
                                                                                   "Gene Set with Specified Rhythmic Pattern" = "rhythmic_gene_set"),
                                                                    inline = TRUE,width = "100%"),
                                                 actionButton(inputId = "retrieve_button",label = "Retrieve"),
                                                 uiOutput(outputId = "bed_files_download")
                                                 
                                                 
                                                 
                                                 
                                                 
                                        )
                            )
                     )
                     
                   )
  )
  
)
