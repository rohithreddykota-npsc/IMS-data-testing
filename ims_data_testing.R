library(moments)
library(ggplot2)
sample_rate <- 20480
nq_freq <- sample_rate/2
rpm = 2000
#rpm = 987
rps = rpm/60
centre_freq <- 236
shaft_freq <- 33.3

imf_path <- "C:/Users/Owner/Documents/ims_data_files"
imf_files <- list.files(imf_path)
imf_full_names <- list.files(imf_path, full.names = T)
raw_data_path <- "D:/sudheer/test2/2nd_test/2nd_test"
raw_files <- list.files(raw_data_path, full.names = T)
axes <- c("x", "y", "z")
fault_ratio <- matrix(0, nrow = 986, ncol = 3)
shaft_ratio <- matrix(0, nrow = 986, ncol = 3)
for (i in 1:length(imf_files)){
  file_name <- imf_files[i]
  file_num <- as.integer(substr(file_name, 10, nchar(file_name)-6))
  axis <- which(axes == substr(file_name, nchar(file_name)-4, nchar(file_name)-4))
  raw_data <- read.table(raw_files[file_num], sep="\t")
  raw_data <- raw_data[,axis]
  raw_data <- raw_data - mean(raw_data)
  imf_data <- read.csv(imf_full_names[i])
  imf_data <- subset(imf_data, select = -X)
  selected_imf <- imf_selection_sum(imf_data, raw_data)
  selected_imf_hilb=abs(hilbert(selected_imf, sample_rate, fftw= FALSE))
  
  fft_data=(2/length(raw_data))*abs(fft(selected_imf))
  fft_data=fft_data[1:floor(length(raw_data)/2)]
  #imf_data_fft <- imf_data_fft[1:2000]
  
  hilb_data=(2/length(raw_data))*abs(fft(selected_imf_hilb))
  hilb_data=hilb_data[1:floor(length(raw_data)/2)]
  
  freqs <- (0:(length(hilb_data)-1))*sample_rate/length(raw_data)
  hilb_data <- as.data.frame(cbind(hilb_data, freqs))
  names(hilb_data) <- c("hilbert", "freqs")
  fault_harmonics <- rbind(c(centre_freq-4, centre_freq + 4),
                           c(centre_freq*2 - 4, centre_freq*2 + 4), 
                           c(centre_freq*3 -4, centre_freq*3 + 4))
  shaft_harmonics <- rbind(c(shaft_freq - 4, shaft_freq + 4), 
                           c(shaft_freq*2 - 4, shaft_freq*2 + 4), 
                           c(shaft_freq*3 - 4, shaft_freq*3 + 4))
  rms_fault <- 0
  rms_shaft <- 0
  for(j in 1:nrow(fault_harmonics)){
    rms_fault <- rms_fault + mean((hilb_data$hilbert[hilb_data$freqs > fault_harmonics[j,1] & 
                                                 hilb_data$freqs < fault_harmonics[j,2]])^2)
    rms_shaft <- rms_shaft + mean((hilb_data$hilbert[hilb_data$freqs > shaft_harmonics[j,1] & 
                                                hilb_data$freqs < shaft_harmonics[j,2]])^2)
  }
  rms_fault <- sqrt(rms_fault/nrow(fault_harmonics))
  rms_shaft <- sqrt(rms_shaft/nrow(fault_harmonics))
  rms_all <- sqrt(mean(hilb_data$hilbert^2))
  fault_ratio[file_num, axis] <- rms_fault/rms_all
  shaft_ratio[file_num, axis] <- rms_shaft/rms_all
  print(i)
}

shaft_x_ratio <- data.frame(cbind(shaft_ratio[1:262,1], 1:262))
names(shaft_x_ratio) <- c("ratio", "file_num")
shaft_x_ratio$type <- "shaft"
fault_x_ratio <- data.frame(cbind(fault_ratio[1:262,1], 1:262))
names(fault_x_ratio) <- c("ratio", "file_num")
fault_x_ratio$type <- "fault"
ratios <- rbind(shaft_x_ratio, fault_x_ratio)
ggplot(ratios, aes(x = file_num, y = ratio, color = type)) + geom_line() + geom_smooth()
ggplot(fault_x_ratio, aes(x = file_num, y = ratio)) + geom_line()

ratio_fft <- abs(fft(shaft_x_ratio$ratio))
ratio_fft <- data.frame(cbind(ratio_fft[2:262], 2:262))
names(ratio_fft) <- c("ratio", "x")
ggplot(ratio_fft, aes(x = x, y = ratio)) + geom_line()

kurtosis(shaft_x_ratio$ratio)
kurtosis(fault_x_ratio$ratio)
