  devtools::document('./Documents/LRAnalytics/')
  devtools::install('./Documents/LRAnalytics/')
  devtools::build('./Documents/LRAnalytics/')


  library(LRAnalytics)
  require(ggplot2)
  paths <- c('SHAM' = '/home/nagai/Documents/Nadine/Report/shan_filtered.csv',
             'IRI' = '/home/nagai/Documents/Nadine/Report/iri_filtered.csv',
             'SNTX' = '/home/nagai/Documents/Nadine/Report/stnx_filtered.csv')


  data <- generate_report(paths,'/home/nagai/Documents/')

