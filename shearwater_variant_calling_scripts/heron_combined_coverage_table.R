
run_2_file_list <- list.files("run_2/")
run_2_file_list <- run_2_file_list[grepl(pattern = "count.csv", x = run_2_file_list)]

run_2_coverage_data_frame <- array(data = NA, dim = c(1181,29903,12))

for(i in 1:1181){
  temp_data <- read.table(paste0("run_2/",run_2_file_list[i]), header = T, sep = ",")
  for(j in 1:12){
    run_2_coverage_data_frame[i,,j] <- temp_data[,j] 
  }
}

 