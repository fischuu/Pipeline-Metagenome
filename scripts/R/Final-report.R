last_n_char <- function(x, n=1){
  substr(x, nchar(x) - n + 1, nchar(x))
}

if(last_n_char(pipelineFolder)!="/") pipelineFolder <- paste0(pipelineFolder, "/")

createRMD.command <- paste0("cat ",pipelineFolder,"scripts/Rmd/final-header.Rmd ",
                                     "> ",projFolder,"/finalReport.Rmd")


cat(createRMD.command, "\n")

system(createRMD.command)

rmarkdown::render(file.path(projFolder,"finalReport.Rmd"), output_file=file.path(projFolder,"finalReport.html"))