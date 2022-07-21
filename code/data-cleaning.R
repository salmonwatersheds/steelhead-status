library(openxlsx)
library("rio")
xls <- dir(pattern = "xlsx")
created <- mapply(convert, xls, gsub("xlsx", "csv", xls))

## Familiarization with dataset

file.info("~/YourDirectoryHere/file_name.csv")$size

#an initial look at the data frame
str(df)


## Check for structural errors

#Variable labels

df <- df %>% rename(employees = How.many.employees.does.your.company.or.organization.have.)

colnames(df)[2]

## Check for data irregularities 

## Document data versions and changes made 