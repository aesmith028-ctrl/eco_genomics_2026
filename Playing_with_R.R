# playing in R to understand basic R functions and objects ####
x <- 5
students <- data.frame(
  name = c("A", "B", "C"),
  height = c(62, 68, 72)
)


head(students)
class(students)
str(students)
# in this list (names), give me the first thing
students$name[1]
# give me the first row, second column
students[1,2]

mean(students$height)
