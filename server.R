# Sys.setlocale("LC_ALL","Russian_Russia.20866")

library(stringr)
library(RDML)
library(ggplot2)
library(tools)
library(plyr)
library(dplyr)
library(MBmca)
library(TSclust)
library(ggdendro)
library(knitr)
library(shiny)

shinyServer(function(input, output, server) {
  values <- reactiveValues(
    hclustResults = list()
  )
  
  rdml <- reactive({
    if (is.null(input$rdmlFile))
      return(NULL)
    isolate({
      newpath <- paste(input$rdmlFile$datapath[1],
                       file_ext(input$rdmlFile$name[1]),
                       sep = ".")
      file.rename(input$rdmlFile$datapath[1],
                  newpath)
      rdml <- RDML$new(newpath)
      description <- rdml$AsTable(name.pattern = paste(react$Position(run$pcrFormat), react$sample$id, data$tar$id))
      mpoints <- rdml$GetFData(description,
                               dp.type = "mdp")
      n.fdata <- ncol(mpoints)
      progress.step <- 1 / n.fdata
      withProgress(
        message = "Сглаживание и нормализаци\u044F образцов",
        value = 0, {
          mpoints.smooth <-    
            cbind(tmp = mpoints[, 1],
                  apply(mpoints[, -1], 2,
                        function(y) {
                          incProgress(progress.step)
                          mcaSmoother(mpoints[, 1],
                                      y,
                                      bgadj = TRUE,
                                      bg = c(20, 30),
                                      minmax = TRUE,
                                      # Trange = trange,
                                      df.fact = 0.6
                          )$y.sp
                          # mcaSmoother(mpoints[, 1],
                          #             y,
                          #             # minmax = TRUE,
                          #             # Trange = trange,
                          #             df.fact = 0.6
                          # )$y.sp
                        })
            ) %>% 
            as.data.frame()
        }
      )
      description$exp.id <- "smooth"
      rdml$SetFData(mpoints.smooth,
                    description, fdata.type = "mdp")
      
      withProgress(
        message = "Вычисление производной",
        value = 0, {
          mpoints.deriv <-    
            cbind(tmp = mpoints.smooth[, 1],
                  apply(mpoints.smooth[, -1], 2,
                        function(y) {
                          incProgress(progress.step)
                          diffQ(
                            data.frame(mpoints.smooth[, 1], y),
                            verbose = TRUE,
                            warn = FALSE)$xy$`d(F) / dT`
                        })
            ) %>% 
            as.data.frame()
        }
      )
      description$exp.id <- "deriv"
      rdml$SetFData(mpoints.deriv,
                    description, fdata.type = "mdp")
      
      mpoints.substr <- as.data.frame(mpoints.deriv)
      progress.step <- 1 / length(unique(description$target))
      withProgress(
        message = "Определение генотипов",
        value = 0, {
          for (test.target in unique(description$target)) {
            fnames.unkn <- (description %>%
                              filter(target == test.target,
                                     sample.type == "unkn"))$fdata.name
            fnames.pos <- (description %>%
                             filter(target == test.target,
                                    sample.type == "pos"))$fdata.name
            deucl <- diss(mpoints.deriv[, c(fnames.unkn, fnames.pos)],
                          "EUCL")
            # zz[[test.target]] <<- deucl 
            values$hclustResults[[test.target]] <- hclust(deucl, "complete")
            deucl <- deucl %>% as.matrix()
            for (fname.unkn in fnames.unkn){
              result <- sapply(fnames.pos, function(fname.pos) {
                deucl[fname.unkn, fname.pos]
                # cor(mpoints.deriv[, fname.unkn],
                #     mpoints.deriv[, fname.pos])
              })
              names(result) <- fnames.pos
              substructed <- mpoints.deriv[, fname.unkn] - mpoints.deriv[, names(which.min(result))]
              
              for (result.name in names(result)) {
                annotation.name <- paste(fname.unkn,  result.name, sep = "_")
                rdml$sample[[description[fname.unkn, "sample"]]]$
                  annotation[[annotation.name]] <- 
                  annotationType$new(
                    annotation.name,
                    as.character(result[result.name])
                  )
#                 values$distances <- rbind(
#                   values$distances,
#                   data.frame(target = test.target,
#                              sample = )
#                 )
              }
              rdml$sample[[description[fname.unkn, "sample"]]]$
                annotation[[paste(fname.unkn,  "distance", sep = "_")]] <- 
                annotationType$new(
                  paste(paste(fname.unkn,  "distance", sep = "_")),
                  as.character(min(result))
                )
              rdml$sample[[description[fname.unkn, "sample"]]]$
                annotation[[paste(fname.unkn,  "fname.pos", sep = "_")]] <- 
                annotationType$new(
                  paste(paste(fname.unkn,  "fname.pos", sep = "_")),
                  names(which.min(result))
                )
              
              mpoints.substr[, fname.unkn] <- substructed
              
            }
          }
          incProgress(progress.step)
        }
      )
      description$exp.id <- "substructed"
      rdml$SetFData(mpoints.substr,
                    description, fdata.type = "mdp")
      rrdml <<- rdml
      rdml
    })
  })
  
  output$targetViewSelector <-  renderUI({
    if (is.null(rdml()))
      return(NULL)
    selectizeInput("viewTarget",
                   "Набор",
                   names(rdml()$target),
                   names(rdml()$target)[1])
  })
  
  output$sampleViewSelector <- renderUI({
    if (is.null(input$viewTarget) ||
        is.null(rdml()))
      return(NULL)
    samples <- (rdml()$AsTable() %>% 
                  filter(target == input$viewTarget,
                         sample.type == "unkn"))$sample %>% 
      unique()
    selectInput("viewSample",
                "Образец",
                samples,
                selected = samples[1]
                , selectize = FALSE
    )
  })
  
  output$hclustPlot <- renderPlot({
    if (is.null(input$viewTarget))
      return(NULL)
    # hcr <<- values$hclustResults[[1]]
    # tt <<- totalTable()
    # source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    
    # colored dendrogram
    # op = par(bg = "#EFEFEF")
    # A2Rplot(values$hclustResults[[input$viewTarget]], k = 3, boxes = FALSE, col.up = "gray50",
    #         col.down = c("#FF6B6B", "#4ECDC4", "#556270"))
    # plot(values$hclustResults[[input$viewTarget]])
    ggdendrogram(values$hclustResults[[input$viewTarget]],
                 rotate = TRUE, color = "red")#, size = 2)
  })
  
  output$meltPlot <- renderPlot({
    if (is.null(input$viewSample) || input$viewSample == "" ||
        is.null(input$viewTarget))
      return(NULL)
    description <- rdml()$AsTable(
      name.pattern = paste(react$Position(run$pcrFormat), react$sample$id, data$tar$id)) %>% 
      filter(target == input$viewTarget)
    
    description.unkn <- description %>% 
      filter(sample == input$viewSample,
             exp.id == "deriv")
    mpoints <-adply(description.unkn$fdata.name, 1, 
                    function(fname) {
                      mpoints.unkn <- rdml()$GetFData(
                        description.unkn %>% 
                          filter(fdata.name == fname),
                        long.table = TRUE,
                        dp.type = "mdp") %>% 
                        mutate(distance = 0)
                      
                      
                      
                      description.pos <- description %>% 
                        filter(exp.id == "deriv",
                               sample.type == "pos")
                      distance.pos <- sapply(description.pos$fdata.name,
                                      function(fname.pos) {
                                        rdml()$sample[[description.unkn[description.unkn$fdata.name == fname, "sample"]]]$
                                          annotation[[paste(str_sub(fname, 1, -2),
                                                            # description.unkn[description.unkn$fdata.name == fname, "target"],
                                                            str_sub(fname.pos, 1, -2), sep = "_")]]$value %>% 
                                          as.numeric()
                                      })
                      description.pos <- description.pos %>% 
                        group_by(fdata.name) %>% 
                        mutate(distance = distance.pos[fdata.name]) %>% 
                        as.data.frame()
                      mpoints.pos <- rdml()$GetFData(
                        description.pos,
                        long.table = TRUE,
                        dp.type = "mdp")
                      mpoints.fname <- bind_rows(mpoints.unkn,
                                                 mpoints.pos)
                      
                      description.substructed <- description %>% 
                        filter(exp.id == "substructed",
                               grepl(str_sub(fname, 1, -2) %>% 
                                       str_replace("\\(", "\\\\("),
                                     fdata.name))
                      # sample == input$viewSample)
                      mpoints.substructed <- rdml()$GetFData(
                        description.substructed,
                        long.table = TRUE,
                        dp.type = "mdp") %>%
                        mutate(sample.type = "substr",
                               distance = 0)
                      mpoints.fname <- 
                        bind_rows(mpoints.fname,
                                  mpoints.substructed) %>%
                        mutate(facet = sprintf("%s    min(r) = %f",
                                               str_sub(fname, 1, -2),
                                               min(distance.pos)))
                      
                      
                      # mpoints.fname <- mpoints.fname %>% 
                      #   group_by(fdata.name) %>% 
                      #   mutate(distance = if (r[1] < 0) distance = 1
                      #          else r)
                      # # rmax <- max(mpoints.fname$r)
                      # # 
                      # # mpoints.fname <- mpoints.fname %>% 
                      # #   filter(sample.type == "unkn") %>% 
                      # #   mutate(distance = rmax)
                      
                      mpoints.fname
                    })
    
    # mpoints <- mpoints %>%
    #   group_by(fdata.name) %>%
    #   mutate(distance = if (distance[1] < 0) distance = 1
    #          else distance)
    # mpoints <- mpoints %>% 
    #   filter(sample.type == "unkn") %>% 
    #   mutate(distance = rmax)
    ggplot(mpoints, 
           aes(x = tmp, y = fluor,
               color = sample.type,
               group = fdata.name,
               alpha = 1 - distance
               # ,linetype = fdata.name
           )) +
      scale_alpha(range=c(0,1), limits=c(0.9, 1), na.value = 0) +
      facet_grid(. ~ facet) +
      geom_line(size = 1.1) +
      geom_point(aes(shape = sample), size = 5)#+
    #coord_cartesian(ylim=c(-1, 1)) 
    # +
    #geom_vline(aes(xintercept = cq))
  })
  
  totalTable <- reactive({
    if (is.null(rdml()))
      return(NULL)
    description <- rdml()$AsTable(
      name.pattern = paste(react$Position(run$pcrFormat), react$sample$id, data$tar$id)) %>% 
      filter(!(exp.id %in% c("smooth", "deriv", "substructed")))
    adply((description %>%
             filter(sample.type == "unkn"))$sample  %>%
            unique, 1,
          function(smpl) {
            adply((description %>%
                     filter(sample.type == "unkn",
                            sample == smpl))$target  %>%
                    unique, 1,
                  function(target) {
                    genotypes <- 
                      sapply(description[description$sample == smpl &
                                           description$target == target &
                                           description$sample.type == "unkn",
                                         "fdata.name"],
                             function(fname) {
                               rdml()$sample[[smpl]]$
                                 annotation[[paste0(fname, "_fname.pos")]]$value
                             }
                      )
                    
                    distances <- 
                      sapply(description[description$sample == smpl &
                                           description$target == target,
                                         "fdata.name"],
                             function(fname) {
                               rdml()$sample[[smpl]]$
                                 annotation[[paste0(fname, "_distance")]]$value %>% 
                                 as.numeric()
                             }
                      )
                    distance <- mean(distances)
                    
                    genotype <- {
                      if (distance > input$distanceLimit) {
                        paste("Генотип не определён, рассто\u044Fние >", input$distanceLimit)
                      } else if (length(unique(genotypes)) > 1 ||
                                 !all(distances <= input$distance.limit)) {
                        "Повторы не совпадают!"
                      } else {
                        description[description$fdata.name == genotypes[1],
                                    "sample"]
                      }
                    }
                    # cat(paste(smpl, target, genotype, distances, sep = "\n"))
                    subresult <- 
                      data.frame(smpl, target, genotype, distance)
                    colnames(subresult) <- c("Образец",
                                          "Набор",
                                          "Генотип",
                                          "Рассто\u044Fние")
                    subresult
                    
                  },
                  .id = "Образец")
          },
          .id = "Образец")
  })
  
  output$resultsTable <- renderDataTable({
    if (is.null(totalTable()))
      return()
    if (input$tabset == "detail")
      totalTable()
    else
      totalTable()[, -4]
  })
  
  output$createReport <- downloadHandler(
    filename  = "report.docx",
    content <- function(file) {
#       src <- normalizePath('report.Rmd')
# #       
# #       # temporarily switch to the temp dir, in case you do not have write
# #       # permission to the current working directory
#       print(tempdir())
#       owd <- setwd(tempdir())
#       print(owd)
#       print(getwd())
#       # on.exit(setwd(owd))
#       file.copy(src, 'report.Rmd')
#       file.copy(src, 'style.docx')
      # knit(input = paste0(getwd(), "/report.Rmd"), 
      #              output = paste0(getwd(), "/report.md"), quiet = FALSE)
      knit(input = "report.Rmd", 
           output = "report.md", quiet = TRUE)
      on.exit(unlink(c("report.md", "figure"), recursive = TRUE))
      # pander::Pandoc.convert(f = paste0(getwd(), "/report.md"), format = "docx",
      #                        options = sprintf("-s --reference-docx=%s/style.docx", getwd()),
      #                        footer = FALSE)
      pander::Pandoc.convert(f =  "report.md", format = "docx",
                             options = "-s --reference-docx=style.docx",
                             footer = FALSE)
#       if (file.exists(paste0(file, ".docx")))
#         file.rename(paste0(file, ".docx"), file)
    }
  )
})
