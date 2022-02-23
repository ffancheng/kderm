# Note: When first deployed, the app will show Error in if: argument is of length zero for the plots location.
# It comes from ml_outlier(x=sim_data()) where the input is reactive
# Just ignore it and wait a few seconds and the plots will show up

library(shiny)
library(shinydashboard)
library(plotly)
library(maniTools)
library(tidyverse)
library(dimRed)
library(viridis)
library(hdrcde)
library(igraph)
library(matrixcalc)
library(ggforce)
library(patchwork)
# library(copula)
Jmisc::sourceAll(here::here("R/sources"))
set.seed(1)

ml_outlier <- function(x, s = 2, k = min(10, nrow(x)), radius = 0, 
                       adjacency = NULL, affinity = NULL,
                       method, 
                       annmethod = c("kdtree", "annoy", "hnsw"),
                       eps = 0, nt = 50, nlinks = 16, ef.construction = 200, ef.search = 10,
                       distance = c("euclidean", "manhattan")[1], diag = FALSE,
                       treetype = c("kd", "bd")[1],
                       searchtype = c("standard", "priority", "radius")[3],
                       perplexity = round(k/3), theta = 0.5, # t-SNE
                       invert.h = TRUE,
                       ell.no = 10, ell.size = 1, gridsize = 20, noutliers = 10,
                       riem.scale = 1, # scaling Riemmanien matrix
                       prob = c(1, 50, 99), # currently not an input of the app
                       ...) {
  
  if (is.null(x) | is.null(method)) return(NULL)
  if(is.null(x$colors)) x$colors <- x$data[,3] # use 3rd column as colors
  
  run_time <- system.time(
    metriclearn <- metricML(x$data, s, k, radius,
                            method = method, 
                            invert.h = invert.h, eps = eps, nt = nt, nlinks = nlinks,
                            annmethod = annmethod, distance = distance, 
                            treetype = treetype,
                            searchtype = searchtype)
  )[[1]]
  
  p_emb <- plot_embedding(metriclearn, color = x$colors, alpha = x$den) +
    labs(x = "", y = "", color = "") +
    plot_ellipse(metriclearn, add = TRUE, ell.no = ell.no,
                 color = blues9[5], fill = blues9[1], alpha = 0.1,
                 ell.size = ell.size) +
    ggtitle(paste0(substring(method, 4), " 2D embedding. Time taken: ", round(run_time, 3), "s.", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))

  fn <- metriclearn$embedding
  Rn <- metriclearn$rmetric
  E1 <- fn[,1]; E2 <- fn[,2]
  # prob <- c(1, 50, 99)
  
  # # TODO: HDR plot for embedding with true density
  # # Convert prob to coverage percentage if necessary
  # if(max(prob) > 50) {# Assume prob is coverage percentage
  #   alpha <- (100 - prob) / 100
  # } else {# prob is tail probability (for backward compatibility)
  #   alpha <- prob}
  # alpha <- sort(alpha)
  # # Calculates falpha needed to compute HDR of bivariate density den.
  # fxy <- x$den
  # falpha <- quantile(fxy, alpha)
  # p <- ggplot2::ggplot(data, ggplot2::aes_string(vnames[1], vnames[2])) +
  #   ggplot2::geom_point(ggplot2::aes_string(col = "Region")) + 
  #   ggplot2::scale_colour_manual(
  #     name = "HDRs",
  #     breaks = c(paste(head(sort(levels), -1)), ">99"),
  #     values = c(RColorBrewer::brewer.pal(length(levels), "YlOrRd")[-1], "#000000")) + 
  #   ggplot2::annotate("text",
  #                      x = data[outliers, 1] + xlim / 50, y = data[outliers, 2] + ylim / 50,
  #                      label = label, col = "blue", cex = 2.5
  # 
  
  f_vkde <- vkde(x = fn, h = Rn*riem.scale, gridsize = gridsize) # estimated densities at grid points
  den <- list(x = f_vkde$eval.points[[1]], y = f_vkde$eval.points[[2]], z = f_vkde$estimate) # VKDE grids
  # fxy <- hdrcde:::interp.2d(f_vkde$x, f_vkde$y, f_vkde$z, x0 = E1, y0 = E2) # linear interpolation for grid estimates
  denhdr <- hdrcde:::den.estimate.2d(x = E1, y = E2, kde.package = "ks", xextend=0.15, yextend = 0.15) # hdr grids
  
  p_vkde <- plot_outlier(x = metriclearn, gridsize = gridsize, prob = prob, noutliers = noutliers, 
    riem.scale = riem.scale, ell.size = ell.size) # f input is for data points, not f_vkde for grids
  p_hdr <- hdrscatterplot_new(E1, E2, kde.package = "ks", noutliers = noutliers)
  p_hdr$p <- p_hdr$p +
    plot_ellipse(metriclearn, ell.no = ell.no, add = TRUE, ell.size = ell.size)

  # p_contour <- filled.contour(f_vkde, color.palette = viridis, # p_vkde$hdr2d_info$den
  #                     plot.axes = { axis(1); axis(2);
  #                       points(fn, pch = 3, col= hcl(c=20, l = 8))})
  # p_contour_hdr <- filled.contour(den, color.palette = viridis, # p_vkde$hdr2d_info$den
  #                     plot.axes = { axis(1); axis(2);
  #                       points(fn, pch = 3, col= hcl(c=20, l = 8))})
  # # p_contour_hdr <- ggplot(as.data.frame(fn), aes(E1, E2)) + geom_density2d_filled()

  
  return(list(metriclearn = metriclearn, p_emb = p_emb, p_vkde = p_vkde, 
              p_hdr = p_hdr, run_time = run_time, estimate = p_vkde$densities, 
              # p_contour = p_contour, p_contour_hdr = p_contour_hdr, 
              den = den, denhdr = denhdr # grid points
              )) #p_contour = p_contour
  
}

# Plotly title font and size
t1 <- list(
  family = "Times New Roman",
  size = 16,
  color = "Black")

algorithm_list <- list("ISOMAP" = "annIsomap",
                       "Locally Linear Embedding (LLE)" = "annLLE",
                       "Laplacian Eigenmaps" = "annLaplacianEigenmaps",
                       # "Hessian LLE" = "annHLLE", 
                       "t-SNE" = "anntSNE",
                       "UMAP" = "annUMAP")

# Define UI for application that plot 2d metadata, 3d dataset, and 2d embedding, HDR plot for outliers
ui <- dashboardPage(

  dashboardHeader(title = "Outlier detection via manifold learning",
                  titleWidth = 350,
                  tags$li(a(href = 'https://github.com/ffancheng/kderm', # private repo
                            icon("github"),
                            title = "Github Repo"),
                          class = "dropdown"
                  )
  ),

  dashboardSidebar(width = 150,
                   sidebarMenu(
                     menuItem("Manifold learning", tabName = "manifoldlearning", icon = icon("table")),
                     menuItem("Comparison", tabName = "comparison", icon = icon("chart-bar"))
                   )
  ),

  dashboardBody(
    tabItems(
      tabItem(
        tabName = "manifoldlearning",
        fluidRow(
          column(4, plotOutput("plot_2d")),
          column(4, plotlyOutput("plot_3d")),
          column(4, plotlyOutput("plot_embed"))
          ),
        br(),

        fluidRow(
          column(8,
            fluidRow(
                   column(6, plotOutput("plot_contour")),
                   column(6, plotlyOutput("plot_vkde"))
                 ),
            br(),
            fluidRow(
              column(6, plotOutput("plot_contour_hdr")),
              column(6, plotlyOutput("plot_outlier"))
              )
          ),
          
          column(4,
                 wellPanel(
                   style = "background-color: #00CC96;",
                   h4("Manifold"),
                   ## Load manifold from files
                   # fluidRow(
                   #   column(4,
                   #          fileInput("file_input", "Load File",
                   #                    accept = c(
                   #                      'text/csv', 'text/comma-separated-values',
                   #                      'text/tab-separated-values', 'text/plain', '.csv', '.tsv')
                   #          ),
                   #          checkboxInput('header', 'Header', TRUE)
                   #   ),
                   #   column(4, radioButtons('sep', 'Separator',
                   #                          c(Comma=',', Semicolon=';', Tab='\t'),
                   #                          ',')
                   #   ),
                   #   column(4, radioButtons('quote', 'Quote',
                   #                          c(None='', 'Double Quote'='"', 'Single Quote'="'"),
                   #                          '"')
                   #   )
                   # ),
                   # # textOutput("file_text"),
                   # hr(),
                   
                   # Generate data from different mappings
                   fluidRow(
                     column(4,
                            selectInput("data_input", "Examples",
                                        choices = c("Swiss Roll", "semi-sphere", "Twin Peak", "S Curve"),
                                        # c("Swiss Roll", "Swiss Hole", "Corner Planes",   
                                        #           "Punctured Sphere", "Twin Peaks", "Clusters",
                                        #           "Toroidal Helix", "Gaussian"),
                                        selected = "Twin Peak")
                     ),
                     column(3, numericInput("num_pts", "No. of Points", 200)),
                     # column(3, uiOutput("data_par")),
                     column(5, selectInput("meta_input", "Meta data type", 
                                           choices = c("copula", "gaussian", "uniform"),
                                           selected = "gaussian"))
                     # column(3, submitButton("Load Example"))  # press button when changing any input
                   ),
                   
                 ),
                 
                 wellPanel(
                   style = "background-color: #636EFA;",
                   h4("Manifold learning parameters"),                            
                   fluidRow(
                     column(3, numericInput("d", "Target Dimension d", 2, min = 2, max = 2)),
                     column(3, selectInput("distance", "Distance", choices = c("euclidean", "manhattan"), selected = "euclidean")),
                     column(3, selectInput("searchtype", "Search Type", choices = list("KNN" = "standard", "Radius" = "radius"), selected = "radius")), 
                     column(3, uiOutput("k_radius")),
                     # column(4, numericInput("k", "Nearest Neighbors k", 8, min = 1)),  # add radius search
                     # column(4, selectInput("kernel", "Kernel",                 
                     #                       choices = c("rbfdot", #"polydot", "vanilladot", "tanhdot",
                     #                                   "laplacedot", "besseldot", "anovadot"#, "splinedot"
                     #                       )))
                   ),
                   # textOutput("comment_text")
                   h5("Approximate nearest neighbor searching"),                                
                   fluidRow(
                     column(6, selectInput("ann", label = "Method", choices = list("k-d trees" = "kdtree", "Annoy" = "annoy", "HNSW" = "hnsw")), seleted = "k-d trees"),
                     column(6, uiOutput("annpar")),   
                   )
                 ),
                 
                 wellPanel(
                   style = "background-color: #4d94ff;",
                   h4("Manifold learning algorithm"),                                  
                   fluidRow(
                     radioButtons("algor", label = "",
                                  choices = algorithm_list,
                                  inline = TRUE)
                   ),
                   # textOutput("plot_text")  ## TODO: not showing the number
                 ),
                 
                 wellPanel(
                   style = "background-color: #ff6666;",
                   h4("Embedding and outlier plot parameters"),
                   fluidRow(
                     column(4, numericInput("ell.no", "Number of ellipses", 10, min = 1, step = 1)),
                     column(4, numericInput("ell.size", "Ellipse size", 0, min = 0, max = 100)),
                     column(4, numericInput("noutliers", "Number of outliers", 10, min = 1, max = 1e8, step = 1)),
                   ),
                   
                   fluidRow(
                     column(4, numericInput("gridsize", "Number of grid points for KDE", 20, min = 1, max = 500, step = 5)),
                     column(4, numericInput("riem.scale", "Riemmanien matrix scaling", 1, min = 0, max = 100)),
                     column(4, checkboxInput("invert.h", label = "Inversed Riemmanien", TRUE))
                   )
                 )
                 
          )       
        )
      ),
      
      
      # Second tab
      tabItem(
        tabName = "comparison",
        fluidRow(
          column(12, offset = 0,
          checkboxGroupInput("algor_group", label = h3("Algorithms"),
                             choices = algorithm_list,
                             inline = TRUE)
        )),
        actionLink("selectall", "Select All"), 
        fluidRow(
          column(1, actionButton("button", "Update")),
          column(11, textOutput("info_text"))
        ),
        br(),
        # plotlyOutput("plot_3d"), # same plot requires different id
        
        fluidRow( uiOutput("plot_outlier_row") )
      )

    )
  )
)




server <- shinyServer(function(input, output, session) {
  
  # # read data from input file and assign colors
  # data_from_file <- reactive({
  #   inFile <- input$file_input
  #   if (is.null(inFile)) return（NULL）
  #   sim_data <- read.csv( inFile$datapath, header = input$header,
  #                         sep = input$sep, quote = input$quote )
  #   if (ncol(sim_data) >= 4) {
  #     colors <- sim_data[,4]
  #   } else {
  #     colors <- z # sim_data[,3], use 3rd dimension as colors  ## TODO, check z
  #   }
  #   list(data = as.matrix(sim_data), colors = colors)
  # })
  
  # reduce_to_3d <- reactive({
  #   sim_data <- data_from_file()
  #   sim_data$data <- sim_data$data[,1:3]
  #   sim_data
  # })

  # DR_data <- reactiveValues(simulation = sim_data())
  # total_time <- reactiveValues(time_taken = NULL)

  output$annpar <- renderUI({
    ann_param_label <- switch(input$ann,
                              "kdtree" = "epsilon",
                              "annoy" = "ntrees",
                              "hnsw" = "nlinks")
    ann_initial_value <- switch(input$ann,
                                "kdtree" = 0,
                                "annoy" = 50,
                                "hnsw" = 16)
    numericInput("ann_parameter", label = ann_param_label, value = ann_initial_value)
  })
  
  output$k_radius <- renderUI({
    searchtype_label <- switch(input$searchtype,
                              "standard" = "No. of nearest neighbors K",
                              "radius" = "Search radius")
    searchtype_initial_value <- switch(input$searchtype,
                                       "KNN" = 10,
                                       "radius" = 10)
    numericInput("searchtype_parameter", label = searchtype_label, value = searchtype_initial_value)
  })
  
  sim_data <- reactive({ mldata(N = input$num_pts, p = 2, meta = input$meta_input, mapping = input$data_input) })
  
  metriclearn_data <- reactive({
    ml_outlier(x = sim_data(), s = input$d, 
               k = input$searchtype_parameter,
               radius = input$searchtype_parameter,
               method = input$algor, 
               invert.h = input$invert.h, eps = input$ann_parameter,
               nt = input$ann_parameter, nlinks = input$ann_parameter,
               annmethod = input$ann, distance = input$distance, 
               treetype = "kd",
               searchtype = input$searchtype,
               ell.no = input$ell.no,
               gridsize = input$gridsize,
               noutliers = input$noutliers, 
               ell.size = input$ell.size,
               riem.scale = input$riem.scale)
  })
  
  # First tab ================================================================================
  output$plot_3d <- renderPlotly({
    # if (!is.null(data_from_file())) {
    #   sim_data <- reduce_to_3d()
    # } else {
      # sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
      # DR_data$simulation <- sim_data
    # }
    if (is.null(sim_data()$data) | (ncol(sim_data()$data) < 3)) {
      plotly_empty(type = "scatter", mode = "markers")
    } else {
      # plotly_3D(sim_data()) %>% layout(title = paste("3D", input$data_input, "mapping data"))
      plot_ly( data.frame(sim_data()$data), x = ~ x, y = ~ y, z = ~ z, color = sim_data()$colors, colors = scales::hue_pal()(4),
               type = "scatter3d", mode = "markers", marker = list(size = 4) ) %>% 
        layout( title = list(text = paste("3D", input$data_input, "mapping data"), font = t1) )
    }
  })
  
  output$plot_2d <- renderPlot({
    # if (!is.null(data_from_file())) {
    #   sim_data <- reduce_to_3d()
    # } else {
      # sim_data <- mldata(N = input$num_pts, meta = input$meta_input, mapping = input$data_input)
      # DR_data$simulation <- sim_data
    # }
    if (is.null(sim_data()$data) | (ncol(sim_data()$data) < 3)) {
      # plotly_empty(type = "scatter", mode = "markers")
      ggplot()
    } else {
      viridis_color <- if(input$meta_input != "gaussian") scale_color_viridis(direction = 1) else NULL # scale_color_manual(values = c('#00CC96', '#EF553B', '#636EFA', '#AB63FA')) # Not matching the default plotly colors
      sim_data()$metadata %>%
        as_tibble() %>% 
        ggplot(aes(x, y, color = sim_data()$colors)) +
        geom_point(aes(alpha = sim_data()$den)) + #alpha = 0.8
        viridis_color + # use four gaussian index for colors
        # scale_color_viridis() +
        labs(color = "Kernels", title = paste("Meta data for", input$data_input, "mapping"), alpha = "Density") +
        theme(plot.title = element_text(hjust = 0.5))
      # plotly_2D(sim_data()$metadata, colors = sim_data()$colors) # if change to plotlyOutput and renderPlotly
    }
  })

  output$plot_embed <- renderPlotly({
    ggplotly(metriclearn_data()$p_emb) %>%
      layout(title = list(font = t1) )
  })

  output$plot_outlier <- renderPlotly({
    ggplotly(metriclearn_data()$p_hdr$p) %>%
      layout(title = list(text = "Outliers with fixed bandwidth", font = t1) )
  })

  output$plot_vkde <- renderPlotly({
    ggplotly(metriclearn_data()$p_vkde$p) %>%
      layout(title = list(text = "Outliers with variable bandwidth", font = t1) )
  })

  output$plot_contour <- renderPlot({
    # metriclearn_data()$p_contour #+
      # ggtitle("Contour plot with variable bandwidth") +
      # theme(plot.title = element_text(family = "Times New Roman", size = 16, color = "Black"))
    filled.contour(metriclearn_data()$den, color.palette = function(n) rocket(n, direction = -1), # p_vkde$hdr2d_info$den
                   plot.title = title(main = "KDE with variable bandwidth", font.main = 1),
                   plot.axes = { axis(1); axis(2);
                     points(metriclearn_data()$metriclearn$embedding, pch = 3, col= sim_data()$colors)})
  })

  output$plot_contour_hdr <- renderPlot({
    filled.contour(metriclearn_data()$denhdr, color.palette = function(n) rocket(n, direction = -1), # p_vkde$hdr2d_info$den
                   plot.title = title(main = "KDE with fixed bandwidth", font.main = 1),
                   plot.axes = { axis(1); axis(2);
                     points(metriclearn_data()$metriclearn$embedding, pch = 3, col= sim_data()$colors)})
  })
  

  # output$file_text <- renderText({"Plot first 3 dimensions only.
  #    Use the 4th dimension as colors if any; otherwise, the 3rd dimension is used."})
  # output$comment_text <- renderText({"The target dimension is fixed at 2."})
  # output$plot_text <- renderPrint({
  #   time <- metriclearn_data()
  #   cat("Time taken:", time$run_time, "s. \n")
  # })
  
  
  # Second tab ================================================================================
  output$info_text <- renderText({"(Note: some algorithms may take long to run (e.g. Isomap and t-SNE),
     please avoid clicking the 'Update' button while the calculation is being performed.)"})

  output$plot_outlier_row <- renderUI({
    plot_output_list <- lapply(1:length(algorithm_list), function(i) {
      fluidRow(
        column(7, plotOutput( paste0("c_plot_", i) )),
        column(5, plotOutput( paste0("f_plot_", i) ), offset = 0)
      )
    })
    do.call(tagList, plot_output_list)
  })


  observeEvent(input$button, {
    algor_list <- input$algor_group
    for (i in 1:length(algorithm_list)) {
      local({
        local_i <- i 
        if ((local_i) %in% seq_along(algor_list)) {
          res <- ml_outlier(x = sim_data(), s = input$d, 
                            k = input$searchtype_parameter,
                            radius = input$searchtype_parameter,
                            method = algor_list[local_i], #input$algor, 
                            invert.h = input$invert.h, eps = input$ann_parameter,
                            nt = input$ann_parameter, nlinks = input$ann_parameter,
                            annmethod = input$ann, distance = input$distance, 
                            treetype = "kd",
                            searchtype = input$searchtype,
                            ell.no = input$ell.no,
                            gridsize = input$gridsize,
                            noutliers = input$noutliers, 
                            ell.size = input$ell.size,
                            riem.scale = input$riem.scale)
        }

        output[[paste0("c_plot_", local_i)]] <-
          renderPlot({
            if ((local_i) %in% seq_along(algor_list)) {
              p3 <- res$p_emb + 
               labs(alpha = "Density", color = "Kernels") + 
               theme(legend.position = "right", legend.direction = "vertical")
              p4 <- res$p_hdr$p + 
                    coord_fixed() + 
                    labs(title = "Fixed bandwidth")
              p5 <- res$p_vkde$p + 
                    coord_fixed() + 
                    labs(title = "Variable bandwidth") 

              p6 <- (p3 + p5 + p4) + 
                plot_layout(guides = 'collect') + 
                plot_annotation(title = substring(algor_list[local_i], 4), 
                                theme = theme(plot.title = element_text(size = 18, face = "bold"))) 
              p6
            } else { ggplot() } # plotly_empty(type = "scatter", mode = "markers")
          })

        output[[paste0("f_plot_", local_i)]] <-
          renderPlot({
            if ((local_i) %in% seq_along(algor_list)) {

              fxy <- sim_data()$den # true meta data density
              fxy_hdr <- res$p_hdr$densities # kde with fixed bandwidth for data points
              fxy_ml <- res$p_vkde$densities # vkde for data points
              
              p1 <- cbind(fxy, fxy_ml) %>% 
                as_tibble() %>% 
                ggplot(aes(fxy, fxy_ml)) + 
                geom_point(col = sim_data()$colors) +
                # scale_x_log10() + 
                # scale_y_log10() + 
                scale_color_manual(values = scales::hue_pal()(4)) + 
                labs(x = "True density", y = "Kernel density estimate",
                     title = paste("Variable bandwidth. Density correlation:", round(cor(fxy, fxy_ml, method = "spearman"), 3) ))
              p2 <- cbind(fxy, fxy_hdr) %>% 
                as_tibble() %>% 
                ggplot(aes(fxy, fxy_hdr)) + 
                geom_point(col = sim_data()$colors) +
                #scale_x_log10() + 
                #scale_y_log10() + 
                scale_color_manual(values = scales::hue_pal()(4)) + 
                labs(x = "True density", y = "Kernel density estimate",
                     title = paste("Fixed bandwidth. Density correlation:", round(cor(fxy, fxy_hdr, method = "spearman"), 3) ))
              p1 + p2

            } else { ggplot() } # plotly_empty(type = "scatter", mode = "markers")
                          
          })

      }) # local
    } # for loop
  }) # observe button

  observe({
    if(input$selectall == 0) return(NULL) 
    else if (input$selectall%%2 == 0)
    {
      updateCheckboxGroupInput(session, "algor_group", "Algorithms", choices = algorithm_list, inline = TRUE)
    } else {
      updateCheckboxGroupInput(session, "algor_group", "Algorithms", choices = algorithm_list, selected = algorithm_list, inline = TRUE)
    }
  })

}) # server


# Run the application
shinyApp(ui = ui, server = server)
