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
  
  # gnomad api query data
  
  rsdat = reactiveValues(rsq = '', ap = 0.650, ep = 0.350)

  observeEvent(input$searchb, {
    rsdat$rsq = trimws(unlist(strsplit(input$rsidsearch, ',')))[1:10]
    rsdat$ap = input$afrprop
    rsdat$ep = input$eurprop
  })
  
  output$propec = renderText({
    shiny::validate(
      need(rsdat$ap <= 1, "African proportion is above 1"),
      need(0 <= rsdat$ap, "African proportion is below 0"),
      need(rsdat$ep <= 1, "European proportion is above 1"),
      need(0 <= rsdat$ep, "European proportion is below 0"),
      need((rsdat$ap + rsdat$ep == 1), "The sum of the proportions is not equal to 1")
    )
  })
  
  output$rsidout = DT::renderDataTable({
    
    shiny::validate(
      need(rsdat$ap <= 1, 'No output, please input a valid rsID list and ancestry proportions'),
      need(0 <= rsdat$ap, 'No output, please input a valid rsID list and ancestry proportions'),
      need(rsdat$ep <= 1, 'No output, please input a valid rsID list and ancestry proportions'),
      need(0 <= rsdat$ep, 'No output, please input a valid rsID list and ancestry proportions'),
      need((rsdat$ap + rsdat$ep == 1), 'No output, please input a valid rsID list and ancestry proportions')
    )
    
    retdat = rsdat$rsq %>%
      rsdfgen(afrp = rsdat$ap, eurp = rsdat$ep) %>%
      as.data.frame() %>%
      select(glink, ref, alt, gnomadafr, updateaf) %>%
      mutate_at(vars(gnomadafr, updateaf), funs(round(., 5))) %>%
      rename('rsID (gnomAD link)' = glink,
             'gnomAD African AF' = gnomadafr,
             'Target Population AF' = updateaf,
             'Ref' = ref,
             'Alt' = alt)
    
    shiny::validate(
      need(dim(retdat)[1] != 0, 'No output, please input a valid rsID list and ancestry proportions')
    )
    
    retdat
    

  }, escape = FALSE,
  options = list(dom = 't',
                 columnDefs = list(list(className = 'dt-center', targets = '_all'))),
  rownames = FALSE)
  
  
  # output$downloadData <- downloadHandler(
  #   
  #   filename = function() {
  #     paste('AncAdjAF', Sys.Date(), '.csv', sep='')
  #   },
  #   content = function(con) {
  #     write.csv(retdat, con)
  #   }
  # )
    
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