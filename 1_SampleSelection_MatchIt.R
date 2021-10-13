library(MatchIt)
library(optmatch)
meta_data <- read.delim('./Final331MHT1DMHT2Dfilter15AttributeData.txt',
                        row.names = 'dbGaP_Subject_ID')
meta_data <- meta_data[ , c(1, 3, 4, 5, 14)]
m.out <- matchit(MHT2D ~ SEX + AGE + RACE, data = meta_data, method = 'optimal', ratio = 2)
summary(m.out)
plot(m.out, type = 'jitter')
plot(m.out, type = 'hist')
match_75_150 <- match.data(m.out)
final_75_150 <- match_75_150[ , 1:5]
final_75_150$SubID = gsub('GTEX[-](+)', '\\1', final_75_150$SUBJID)
write.table(final_75_150, file = './final_75_150.txt', 
            row.names = FALSE, sep = '\t')