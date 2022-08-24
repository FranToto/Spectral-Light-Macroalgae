##############################################################################
# Shiny app to explore Synth_vars_RGB_v09
#
# FT - 24/08/2022
# Make it tidier for release on GitHub/Shinyapps.io
##############################################################################

library(shiny)
library(shinybusy)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(ggtern)
library(scico)
library(colorspace)

#### Data ####
vars_all <- readRDS('Synth_vars_RGB_v09.rds') %>% 
  mutate(bbp_rank = rank(bbp),
         adet_rank = rank(adet),
         aph_rank = rank(aph),
         bbp_trans = 1 - log(1000000 - bbp_rank +1)/log(1000000),
         adet_trans = 1 - log(1000000 - adet_rank +1)/log(1000000),
         aph_trans = 1 - log(1000000 - aph_rank +1)/log(1000000),
         PercentLightAction_brown = EbedAction_brown*100/Ebed_par,
         PercentLightAction_red_scaled = EbedAction_red_scaled*100/Ebed_par,
         PercentLightAction_green_scaled = EbedAction_green_scaled*100/Ebed_par,
         Ebed_ratio_brown = EbedAction_brown/EbedAction_red_scaled,
         Ebed_ratio_red_scaled = EbedAction_green_scaled/EbedAction_red_scaled,
         Ebed_ratio_green_scaled = EbedAction_green_scaled/EbedAction_brown) 


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Explore Synthetic Dataset v09"),

      add_busy_spinner(timeout=1000,color='blue'),
      splitLayout(style = "border: 1px solid silver:", cellWidths = c("33%", "33%","33%"), 
                plotOutput("brown"),
                plotOutput("red"),
                plotOutput("green")
      ),
    hr(),
    
    fluidRow(
      column(3,
             h4("Enter values"),
             numericInput("z",
                        "Depth:",
                         min = 1,
                         max = 50,
                         step = 0.1,
                         value = 10),
             br(),
             numericInput("kpar",
                          "Kd(PAR):",
                          min = 0.05,
                          max = 1,
                          step = 0.01,
                          value = 0.2),
             
      ),
      column(4, offset = 1,
             radioButtons('metrics', 
                          'Metrics', 
                          choices = c('Performance'='EbedAction',
                                      'Efficiency'='PercentLightAction',
                                      'Competitive Performance'='Ebed_ratio',
                                      'Taxon Dominance'='Taxa_winner'),
                          selected = 'PercentLightAction')
      ),
      column(4,
             numericInput("ptsize",
                         "Size of points:",
                         min = 0.7,
                         max = 3,
                         step = 0.1,
                         value = 2)
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


        output$brown <- renderPlot({
            
            # Add condition to show 3 plots if metric is not Dominance
            if(!input$metrics=='Taxa_winner'){
                
                p <-   vars_all %>%
                    filter(between(depth,input$z - 0.5,input$z + 0.5),between(bbp,0,0.1),between(aph,0,0.2),between(Kpar,input$kpar - 0.01,input$kpar + 0.01)) %>% 
                    ggtern(data=.,aes(x=bbp_trans,y=adet_trans, z=aph_trans)) +
                    geom_point(aes_string(col=paste0(input$metrics,'_brown')),size=input$ptsize) +
                    xlab(expression(b[bp])) + ylab(expression(a[ph])) + zlab(expression(a[det])) +
                    Tarrowlab('aph') + Larrowlab('bbp') + Rarrowlab('adet') +
                    ggtitle(paste0(input$metrics,' Brown ')) +
                    theme_rgbw() +
                    theme(legend.title=element_blank())
                
                # Add condition to change color scale depending on metrics
                if(input$metrics=='Ebed_ratio'){
                    p <- p + scale_colour_continuous_diverging(palette = "Purple-Brown",mid=1) + theme_rgbg(base_size = 12, base_family = "") + theme(legend.title=element_blank())
                }else{
                    p <- p + scale_color_viridis()
                } 
                
                print(p)#Key, thanks to https://stackoverflow.com/questions/55101305/ternary-plot-of-ggtern-not-working-in-shiny
            
            }else{#Only 1 plot for Dominance metrics
                p_dominance <- vars_all %>% 
                    filter(.,between(depth,input$z - 0.5,input$z + 0.5),between(Kpar,input$kpar - 0.01,input$kpar + 0.01)) %>% 
                    ggtern(data=.,aes(x=bbp_trans,y=adet_trans, z=aph_trans)) +
                    geom_point(aes(col=Taxa_winner),size=input$ptsize) +
                    scale_colour_manual(values = c("#DDAA33", "#009988", "#CC3311"),labels = c("Brown", "Green", "Red")) +
                    xlab(expression(b[bp])) + ylab(expression(a[ph])) + zlab(expression(a[det])) +
                    Tarrowlab('aph') + Larrowlab('bbp') + Rarrowlab('adet') +
                    theme_rgbw()+
                    theme(legend.title=element_blank())
                print(p_dominance)
                
            }     
        })
        
        
        output$red <- renderPlot({
            
            # Add condition to show 3 plots if metric is not Dominance
            if(!input$metrics=='Taxa_winner'){
                
                p <- vars_all %>%
                    filter(between(depth,input$z - 0.5,input$z + 0.5),between(bbp,0,0.1),between(aph,0,0.2),between(Kpar,input$kpar - 0.01,input$kpar + 0.01)) %>% 
                    ggtern(data=.,aes(x=bbp_trans,y=adet_trans, z=aph_trans)) +
                    geom_point(aes_string(col=paste0(input$metrics,'_red_scaled')),size=input$ptsize) +#Thanks to https://stackoverflow.com/questions/35345782/shiny-passing-inputvar-to-aes-in-ggplot2
                    xlab(expression(b[bp])) + ylab(expression(a[ph])) + zlab(expression(a[det])) +
                    Tarrowlab('aph') + Larrowlab('bbp') + Rarrowlab('adet') +
                    ggtitle(paste0(input$metrics,' Red ')) +
                    theme_rgbw() +
                    theme(legend.title=element_blank())
                
                if(input$metrics=='Ebed_ratio'){
                    p <- p + scale_colour_continuous_diverging(palette = "Red-Green",mid=1) + theme_rgbg(base_size = 12, base_family = "") + theme(legend.title=element_blank())
                }else{
                    p <- p + scale_color_viridis()
                } 
                print(p)#Key, thanks to https://stackoverflow.com/questions/55101305/ternary-plot-of-ggtern-not-working-in-shiny
                
            }else{#Only 1 plot for Dominance metrics
                print('')
            } 
            
        })
        
        
        output$green <- renderPlot({
            
            # Add condition to show 3 plots if metric is not Dominance
            if(!input$metrics=='Taxa_winner'){
                
                p <- vars_all %>%
                    filter(between(depth,input$z - 0.5,input$z + 0.5),between(bbp,0,0.1),between(aph,0,0.2),between(Kpar,input$kpar - 0.01,input$kpar + 0.01)) %>% 
                    ggtern(data=.,aes(x=bbp_trans,y=adet_trans, z=aph_trans)) +
                    geom_point(aes_string(col=paste0(input$metrics,'_green_scaled')),size=input$ptsize) +
                    xlab(expression(b[bp])) + ylab(expression(a[ph])) + zlab(expression(a[det])) +
                    Tarrowlab('aph') + Larrowlab('bbp') + Rarrowlab('adet') +
                    ggtitle(paste0(input$metrics,' Green ')) +
                    theme_rgbw() +
                    theme(legend.title=element_blank())
                
                if(input$metrics=='Ebed_ratio'){
                    p <- p + scale_colour_continuous_diverging(palette = "Green-Brown",mid=1,rev=T) + theme_rgbg(base_size = 12, base_family = "") + theme(legend.title=element_blank())
                }else{
                    p <- p + scale_color_viridis()
                } 
                print(p)#Key, thanks to https://stackoverflow.com/questions/55101305/ternary-plot-of-ggtern-not-working-in-shiny
                
            }else{#Only 1 plot for Dominance metrics
                print('')
            } 
            
        })
        

}

# Run the application 
shinyApp(ui = ui, server = server)
