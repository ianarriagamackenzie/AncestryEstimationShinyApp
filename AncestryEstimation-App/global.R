# Hendricks Mixtures Research Group - UC Denver
# ver 3
# Created by Ian Arriaga MacKenzie
# For use with gnomAD ancestry analysis information

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(dplyr)
library(shinycssloaders)
library(skimr)

# sets the colors used for the 5 ancestries, consistent coloring across different plots
colvec = c(brewer.pal(n = 8, name = 'Set2')[3],
           brewer.pal(n = 8, name = 'Set2')[1],
           brewer.pal(n = 8, name = 'Set2')[2],
           brewer.pal(n = 8, name = 'Set2')[6],
           brewer.pal(n = 8, name = 'Set2')[5])

# load genome/exome test data
dat = read.csv('findat.txt', sep='')

# load data from laptop for testing purposes, comment out on published app
# dat = read.csv('~/GitHub/mixturesresearch/AncEstTestApp/findat.txt', sep='')

# Functions

# main proportion plot of ancestries for block bootstrap and random snp sample options
proportionplot = function(dataset){
  
  # select correct columns, melt for plotting
  meltdat = dataset %>%
    select(AFR, EAS, EUR, NAM, SAS) %>% 
    rename(IAM = NAM) %>% 
    melt()
  names(meltdat) = c('Ancestry', 'Proportion')
  
  # ggplot with options
  propplot = ggplot(meltdat, aes(x=Proportion, fill = Ancestry)) +
    geom_histogram(bins = 400)+
    facet_grid(Ancestry ~ .)+
    scale_fill_manual(values = colvec)+
    guides(fill = FALSE)+
    scale_x_continuous(breaks = c(0,.25,.50,.75,1),
                       limits = c(-0.01,1)) +
    theme(axis.text.x = element_text(size=20),
          axis.title.x = element_text(size=20),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y = element_text(size = 18),
          axis.ticks = element_blank())
  
  return(propplot)
}

# secondary plot of individual ancestry distributions and 95% CIs, for block bootstrapping and random snp sample
distplot = function(dataset){
  
  # select correct columns, set up return plot list
  datasetc = dataset %>% 
    select(AFR, EAS, EUR, NAM, SAS) %>% 
    rename(IAM = NAM)
  plot_list = list()
  anc_list = names(datasetc)
  
  # ggplot with options, returns list of 5 ancestry plots
  for (i in 1:5){
    pl = ggplot(data = datasetc, aes_string(x = anc_list[[i]])) +
      geom_histogram(aes(y = ..density..), color = 'black', fill = colvec[i], bins = 30) +
      geom_density(fill = NA) +
      geom_vline(data = dataset, xintercept = quantile(datasetc[[anc_list[i]]], probs = 0.025),
                 linetype ="longdash", size = .8) +
      geom_vline(data = dataset, xintercept = quantile(datasetc[[anc_list[i]]], probs = 0.975),
                 linetype ="longdash", size = .8) +
      theme_minimal() + 
      labs(title = paste(anc_list[i]), x = NULL, y = NULL) +
      scale_x_continuous(labels = percent_format(accuracy = .01), 
                         breaks = c(quantile(datasetc[[anc_list[i]]], probs = 0.025),
                                    quantile(datasetc[[anc_list[i]]], probs = 0.975))) +
      theme(plot.title = element_text(hjust = 0.5, size = 15),
            axis.text.x = element_text(face="bold", size=15),
            axis.text.y = element_blank())
    plot_list[[i]] = pl
  }
  
  # return using gridExtra library
  return(grid.arrange(plot_list[[1]],
                      plot_list[[2]],
                      plot_list[[3]],
                      plot_list[[4]],
                      plot_list[[5]], ncol = 5))
}

# plotting ancestry proportions by chromosome
chrplot = function(dataset){
  
  # subset correct columns, name and melt them
  chrdat = dataset %>%
    select(TestType, AFR, EAS, EUR, NAM, SAS)
  names(chrdat) = c('Chromosome', 'AFR', 'EAS', 'EUR', 'IAM', 'SAS')
  chrmelt<- melt(chrdat, id="Chromosome", 
                 measure=c('AFR', 'EAS', 'EUR', 'IAM', 'SAS'), 
                 variable.name="Anc", value.name="Proportions")
  
  # ggplot with options
  chrplot = ggplot(chrmelt, aes(Chromosome, Proportions, fill=Proportions)) + 
    facet_wrap( ~ Anc , nrow = 1) +
    geom_bar(stat="identity") +
    ylim(c(0,1)) +
    coord_flip() +
    scale_fill_distiller(palette = 'Spectral')+
    guides(fill = FALSE) +
    scale_x_reverse(breaks = c(1:22), expand = c(0,0)) +
    theme(
      panel.grid.minor.y = element_blank()
    )
  
  return(chrplot)
  
}
# numeric information plots for block bootstrapping and random snp sample
bbraninfo = function(dataset){
  
  # select and summarize correct columns, return
  # mean, standard deviation, min, 25th percentile, median, 75th percentile, max
  sdat = dataset %>% 
    select(AFR, EAS, EUR, NAM, SAS) %>% 
    rename(IAM = NAM) %>% 
    skim() %>% 
    select(skim_variable, numeric.mean, numeric.sd, numeric.p0, 
           numeric.p25, numeric.p50, numeric.p75, numeric.p100)
  names(sdat) = c('Ancestry', 'Mean', 'SD', 'P0', 'P25', 'P50', 'P75', 'P100')
  return(sdat)
  
}