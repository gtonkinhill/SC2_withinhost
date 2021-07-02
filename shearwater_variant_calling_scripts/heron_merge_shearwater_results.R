files <- list.files("background_panel_2_run_1/")
for(i in 1:length(files)){
  temp_data <- readRDS(paste0("background_panel_2_run_1/",files[i]))
  if(i == 1){
    combined_data <- temp_data
  } else{
    combined_data <- rbind(combined_data,temp_data)
  }
}
saveRDS(object = combined_data, file = "initial_shearwater_calls_background_panel_2_run_1.rds")
