# shiny server, ancestry mixtures research project
server = function(input, output, session) {
  
  # store correct data frames in reactive value
  bbranchrdat = reactiveValues(bbdat = NULL, randat = NULL, chrdat = NULL)
  
  # updates reactive data frame when new gnomAD group or genome/exome is selected
  observeEvent({
    input$ancdat
    input$exge}, {
      
      tempdat = dat %>% 
        filter(Exge %in% paste(input$exge)) %>% 
        filter(Gnomadanc %in% paste(tolower(input$ancdat)))
      
      bbranchrdat$bbdat = tempdat %>% 
        filter(TestType %in% 'blockbootstrap')
      
      bbranchrdat$randat = tempdat %>% 
        filter(TestType %in% 'Random_SNPs_Genomewide')
      
      bbranchrdat$chrdat = tempdat %>% 
        filter(TestType %in% levels(dat$TestType)[grep('CHR', levels(dat$TestType))])
      bbranchrdat$chrdat$TestType = c(1:22)
      
    })
  
  # output block bootstrapping main proportion plot
  output$plotbb = renderPlot({
    
    bbranchrdat$bbdat %>% 
      proportionplot()
    
  })
  
  # output block bootstrapping secondary distribution plot
  output$distbb = renderPlot({
    
    bbranchrdat$bbdat %>% 
      distplot()
    
  })
  
  # output block bootstrapping numeric summary
  output$infobb = renderTable({
    
    bbranchrdat$bbdat %>% 
      bbraninfo()
    
  }, digits = 5, align = 'c')
  
  # output random snp sample main proportion plot
  output$plotran = renderPlot({
    snpnum = 0
    if (input$exge == 'genome'){
      snpnum = as.numeric(input$randsnpnumge)
    }
    else if (input$exge == 'exome'){
      snpnum = as.numeric(input$randsnpnumex)
    }
    
    bbranchrdat$randat %>%
      filter(NumberSNPs == snpnum) %>% 
      proportionplot()
    
  })
  
  # output random snp sample secondary distribution plot
  output$distran = renderPlot({
    snpnum = 0
    if (input$exge == 'genome'){
      snpnum = as.numeric(input$randsnpnumge)
    }
    else if (input$exge == 'exome'){
      snpnum = as.numeric(input$randsnpnumex)
    }
    
    bbranchrdat$randat %>%
      filter(NumberSNPs == snpnum) %>% 
      distplot()
    
  })
  
  # output random snp sample numeric summary
  output$inforan = renderTable({
    snpnum = 0
    if (input$exge == 'genome'){
      snpnum = as.numeric(input$randsnpnumge)
    }
    else if (input$exge == 'exome'){
      snpnum = as.numeric(input$randsnpnumex)
    }
    
    bbranchrdat$randat %>%
      filter(NumberSNPs == snpnum) %>% 
      bbraninfo()
    
  }, digits = 5, align = 'c')
  
  # output chromosome main plot
  output$plotchr = renderPlot({
    
    bbranchrdat$chrdat %>% 
      chrplot()
    
  })
  
  # output chromosome numeric summary
  output$sumchr = renderTable({
    
    sumchrdat = bbranchrdat$chrdat %>% 
      select(TestType, AFR, EAS, EUR, NAM, SAS)
    names(sumchrdat) = c('Chromosome', 'AFR', 'EAS', 'EUR', 'IAM', "SAS")
    
    return(sumchrdat)
    
  }, digits = 5, spacing = 'xs', align = 'c')
  
}