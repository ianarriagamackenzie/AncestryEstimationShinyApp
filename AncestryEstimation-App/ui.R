# shiny ui, ancestry mixtures research project
ui = dashboardPage(
  
  skin = 'blue',
  
  # main header
  dashboardHeader(
    
    title = 'gnomAD Ancestry Estimation (BETA)',
    titleWidth = 365
    
  ),
  
  # options for choices along side bar
  dashboardSidebar(
    width = 250,
    
    # secondary header
    fluidRow(
      align = 'center',
      h4('Hendricks Research Group'),
      h4('University of Colorado Denver')
    ),
    
    # user selection of the exome or genome
    radioGroupButtons(
      inputId = 'exge',
      label = NULL,
      choiceNames = c('Genome', 'Exome'),
      choiceValues = c('genome', 'exome'),
      selected = 'genome',
      individual = TRUE,
      width = '100%',
      justified = TRUE,
      status = 'primary',
      checkIcon = list(
        yes = icon("ok", 
                   lib = "glyphicon"))
    ),
    
    # user selection of gnomAD group
    pickerInput(
      inputId = 'ancdat',
      label = 'Ancestry Group',
      choices = c('AFR', 'AMR', 'OTH'),
      choicesOpt = list(
        subtext = c('African/African American', 'American/Latinx', 'Other')
      ),
      selected = 'AFR',
      options = list(
        `live-search` = TRUE)
    ),
    
    # main menu of different plots
    sidebarMenu(
      id = 'menuselect',
      
      menuItem("Ancestry Proportions",
               icon = icon("chart-area"),
               startExpanded = TRUE,
               
               # Block Boostrap analysis
               menuSubItem("Block Bootstrap",
                           tabName = "bb",
                           selected = TRUE),
               
               # Random SNP sample analysis
               menuSubItem("Random SNP Sample",
                           tabName = "ran")
      ),
      
      # analysis by chromosome
      # menuItem("Ancestry Proportions by Chromosome", tabName = "chr", icon = icon("chart-bar")),
      
      menuItem("Ancestry Adjusted Allele Frequency", tabName = "adjaf", icon = icon("calculator")),
      
      # readme with disclaimer, acknowledgements, general use information
      menuItem("ReadMe", tabName = "readme", icon = icon("readme"))
      
      # github link
      # menuItem("Github", icon = icon("code"),
      #          href = "https://github.com/hendriau/Mixtures",
      #          newtab = TRUE)
      
    )
    
  ),
  
  dashboardBody(
    
    # change the font size to 12
    tags$head(
      tags$style(HTML(".main-sidebar { font-size: 12px; }"))
    ),
    
    tabItems(
      
      # block bootstrapping page
      tabItem(tabName = "bb",
              
              fluidRow(
                
                column(
                  width = 3,
                  
                  # description of block bootstrapping
                  box(
                    title = "Block Bootstrapping", width = NULL, status = "primary",
                    'We use block bootstrapping to estimate error for the ancestry proportions. 
                    We resample 3,357 centiMorgan blocks 1,000 times for the plots and confidence intervals shown here.'
                  )
                  ),
                
                column(
                  width = 9,
                  
                  # main proportion plot panel
                  tabBox(
                    title = "Proportion Estimates for Block Bootstrapping",
                    width = NULL, height = 440, side = 'right', selected = 'Visual',
                    
                    # numeric summary
                    tabPanel(
                      'Numeric',
                      withSpinner(tableOutput(
                        'infobb'
                      ))
                    ),
                    
                    # visual plots
                    tabPanel(
                      'Visual',
                      withSpinner(plotOutput(
                        'plotbb',
                        height = 370
                      ))
                    )
                  )
                )
                
                ),
              
              fluidRow(
                
                column(width = 1),
                column(
                  width = 10,
                  
                  # distribution plot panel
                  box(
                    title = "Distribution Plots and 95% Confidence Intervals", width = NULL,
                    status = "primary", height = 240,
                    withSpinner(plotOutput(
                      'distbb',
                      height = 170
                    ))
                  )
                )
                
              )
      ),
      
      # random SNP sample page
      tabItem(tabName = "ran",
              
              fluidRow(
                
                column(
                  width = 3,
                  
                  # description of random SNP sample analysis
                  box(
                    title = "Random SNP Sample", width = NULL, status = "primary",
                    'We sample N random SNPs across the 22 autosomes to estimate ancestry proportions. 
                    We randomly sample 1,000 times for the plots and confidence intervals shown here. 
                    N can be varied to evaluate our method with different numbers of SNPs.'
                  ),
                  
                  # panel to select number of SNPs randomly sampled
                  box(
                    title = 'N Random SNPs', width = NULL, status = "primary",
                    
                    # genome number selection
                    conditionalPanel(
                      condition = "input.exge == 'genome' ",
                      sliderTextInput(
                        inputId = 'randsnpnumge',
                        label = NULL,
                        choices = c(10, 50, 100, 500, 1000, 2500, 5000, 10000, 50000, 100000),
                        selected = '1000',
                        grid = TRUE,
                        hide_min_max = TRUE
                      )
                    ),
                    
                    # exome number selection
                    conditionalPanel(
                      condition = "input.exge == 'exome' ",
                      sliderTextInput(
                        inputId = 'randsnpnumex',
                        label = NULL,
                        choices = c(10, 50, 100, 500, 1000, 2500, 5000, 7500, 9000),
                        selected = '1000',
                        grid = TRUE,
                        hide_min_max = TRUE
                      )
                    )
                  )
                  
                  ),
                
                column(
                  width = 9,
                  
                  # main proportion plot panel
                  tabBox(
                    title = "Proportion Estimates for Random SNP Sample",
                    width = NULL, height = 440, side = 'right', selected = 'Visual',
                    
                    # numeric summary
                    tabPanel(
                      'Numeric',
                      withSpinner(tableOutput(
                        'inforan'
                      ))
                    ),
                    
                    # visual plots
                    tabPanel(
                      'Visual',
                      withSpinner(plotOutput(
                        'plotran',
                        height = 370
                      ))
                    )
                  )
                )
                
                ),
              
              fluidRow(
                
                column(width = 1),
                column(
                  width = 10,
                  
                  # distribution plot panel
                  box(
                    title = "Distribution Plots and 95% Confidence Intervals", width = NULL,
                    status = "primary", height = 240,
                    withSpinner(plotOutput(
                      'distran',
                      height = 170
                    ))
                  )
                )
                
              )
  ),
  
  # chromosome analysis page
  tabItem(tabName = "chr",
          
          fluidRow(
            
            column(
              width = 3,
              
              # chromosme analysis information panel
              box(
                title = "Chromosome", width = NULL, status = "primary",
                'Estimated ancestry proportions by chromosome using all SNPs.'
              )
            ),
            
            column(
              width = 9,
              
              # main chromosome plot panel
              tabBox(
                title = "Proportion Estimates by Chromosome",
                width = NULL, height = 700, side = 'right', selected = 'Visual',
                
                # numeric summary
                tabPanel(
                  'Numeric',
                  withSpinner(tableOutput(
                    'sumchr'
                  ))
                ),
                
                # visual plot
                tabPanel(
                  'Visual',
                  withSpinner(plotOutput(
                    'plotchr',
                    height = 650
                  ))
                )
              )
            )
            
          )
          
  ),
  
  #ancestry adjusted function tab
  tabItem(tabName = 'adjaf',
          fluidRow(
            
            column(
              width = 4,
              
              box(title = 'Ancestry Adjusted Allele Frequency', width = NULL, status = 'primary',
                  'Using the ancestry proportions estimated in the gnomAD African/African American group, 
                  we adjust the allele frequencies to match a user given target population at specific SNPs',
                  br(),
                  br(),
                  'BETA: gnomAD v2.1.1. AFR group only at present.'
                  ),
              
              box(title = "User Input", width = NULL, status = "primary",
                  
                  'List of rsIDs, comma delimited.',
                  
                  br(),
                  br(),
                  
                  textAreaInput(
                    inputId = "rsidsearch", label = NULL,
                    placeholder = "rs4244285, rs3093024, rs744373, rs597668, ...",
                    rows = 3,
                    resize = 'none',
                    width = "100%"
                  ),
                  
                  'Ancestry proportions of a target admixed population.',
                  
                  br(),
                  br(),
                  
                  splitLayout(
                    numericInput("afrprop", 
                                 'African', 
                                 .650, 
                                 min = 0, max = 1, step = .001),
                    
                    numericInput("eurprop",
                                 'Non-Finnish European',
                                 .350,
                                 min = 0, max = 1, step = .001)
                  ),
                  
                  textOutput('propec'),
                  
                  br(),
                  
                  actionBttn(
                    inputId = "searchb",
                    label = 'Search',
                    style = 'jelly',
                    color = "primary",
                    icon = icon("search")
                  ),
                  
                  br(),
                  br(),
                  
                  'Querying gnomAD database may take several seconds.'
                  
                  
              )
            ),
            
            column(
              width = 8,
              
              box(title = 'Allele Frequency Output Table', width = NULL, status = 'primary',
              
              # downloadLink('downloadData', 'Save Results'),
              
              withSpinner(DT::dataTableOutput('rsidout')),
              tags$style(type="text/css", '#rsidout tfoot {display:none;}')
              
              )
              
            )
          
          )
  ),
  
  # readme page
  tabItem(tabName = "readme",
          fluidRow(
            
            column(
              width = 8,
              
              # main readme, describes data sets and overall process of project, general information
              box(
                title = "ReadMe", width = NULL, status = "primary",
                
                strong('Purpose'),
                br(),
                'Estimate the proportion of reference ancestry groups in summary genotype frequency data.',
                br(),
                br(),
                strong('Data'),
                br(),
                'Our reference panel was created from ',
                a('1000 Genomes Project', href = "https://www.internationalgenome.org/", target="_blank"),
                ' (GRCh37/hg19) superpopulations (African, Non-Finish European, East Asian, South Asian) and an ',
                a('Indigenous American population', href = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130711_native_american_admix_train", target="_blank"),
                ' (616,568 SNPs and 43 individuals, GRCh37/hg19). Tri-allelic SNPs and SNPs with missing 
                allele frequency information were removed, leaving 613,298 SNPs across the 22 autosomes.',
                br(),
                br(),
                'We estimate the ancestry proportions from ',
                a('gnomAD V2', href = 'https://gnomad.broadinstitute.org/', target="_blank"),
                '(GRCh37/hg19). After merging with our reference panel we checked for allele matching and strand flips. 
                Our final dataset had 582,550 genome SNPs and 9,835 exome SNPs across the 22 autosomes.'
              ),
              
              # disclaimer about use of method
              box(
                title = "Disclaimer", width = NULL, status = "primary",
                
                'Under no circumstances shall authors of this website and ancestry estimation algorithm be liable for 
                    any indirect, incidental, consequential, special or exemplary damages arising out of or in connection
                    with your access or use of or inability to access the ancestry estimation website or any associated software 
                    and tools and any third party content and services, whether or not the damages were foreseeable and whether or 
                    not the authors were advised of the possibility of such damages.  By using the ancestry estimation platform 
                    you agree to use it to promote scientific research, learning or health.'
              )
              ),
            
            column(
              width = 4,
              
              # acknowledgements and thank yous, contact information
              box(
                title = "Acknowledgements", width = NULL, status = "primary",
                
                strong("This work was a collaborative effort by:"),
                br(), 
                "Ian S. Arriaga Mackenzie, Gregory M. Matesi, Alexandria Ronco, Ryan Scherenberg, 
                Andrew Zerwick, Yinfei Wu, James Vance, Sam Chen, Kaichao Chang, Katie Marker, Jordan R. Hall, 
                Christopher R. Gignoux, Megan Null, Audrey E. Hendricks",
                br(),
                br(),
                strong("Additional Funding"),
                br(),
                "CU Denver Undergraduate Research Opportunity Program (UROP)",
                br(),
                "CU Denver Education through Undergraduate Research and Creative Activities program (EUReCA)",
                br(),
                br(),
                # shiny app creater, Ian Arriaga MacKenzie, contact information
                strong('Shiny App'),
                br(),
                'Ian S. Arriaga MacKenzie',
                br(),
                'IAN.ARRIAGAMACKENZIE@ucdenver.edu',
                br(),
                a(actionButton(inputId = "email1", label = NULL,
                               icon = icon("envelope", lib = "font-awesome")),
                  href="mailto:IAN.ARRIAGAMACKENZIE@ucdenver.edu"),
                br(),
                br(),
                # PI for project, Dr. Audrey Hendricks, contact information
                strong('Principal Investigator'),
                br(),
                'Audrey E. Hendricks, Ph.D.',
                br(),
                'AUDREY.HENDRICKS@ucdenver.edu',
                br(),
                a(actionButton(inputId = "email2", label = NULL,
                               icon = icon("envelope", lib = "font-awesome")),
                  href="mailto:AUDREY.HENDRICKS@ucdenver.edu")
              ),
              
              # CU Denver logo
              img(src='CUdenverlogo.png', align = "Center", height = 150, width = 270)
            )
            
            )
          
  )
  
)

)

)