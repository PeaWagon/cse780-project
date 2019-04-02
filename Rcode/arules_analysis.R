
# association rules analysis

library(arules)
source("std_lift.R") # professor's code to calculate
                     # standardised lift

infile <- "../data_processing/drug_data_trimmed_text_len.csv"

d <- read.csv(infile)

# we can consider:
# drugName (1), conditon (2), rating (4), revLenDesc(9)
x <- d[, c(1,2,4,9)]

# put the rating in two categories: >= 7 and < 7
# sets rating variable to 0 for values below 7
x$rating[x$rating<7] <- 0

# sets rating variable to 1 for values above or equal to 7
x$rating[x$rating>0] <- 1

# make the option 0 or 1
x$rating <- as.factor(x$rating)

# minlen is number of parameters to consider
# maxlen is how many parameters there are for the rules
# if no rules are generated, reduce support and/or
# confidence
params <- list(support=0.01, confidence=0.8, minlen=2,
               maxlen=4)

# set the right-hand side or the left-hand side
# to see what rules there are (a)
app <- list(rhs=c("rating=0", "rating=1"),
            default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# see where drugs are poorly rated (b)
params <- list(support=0.005, confidence=0.6, minlen=2,
               maxlen=4)

app <- list(rhs=c("rating=0"), default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# see if length of text is a rule for anything (c)
params <- list(support=0.0001, confidence=0.7,
               minlen=2, maxlen=4)

app <- list(rhs=c("revLenDesc=v.long",
                  "revLenDesc=long",
                  "revLenDesc=medium",
                  "revLenDesc=short",
                  "revLenDesc=v.short"), default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# remove rhs restrictions (d)
params <- list(support=0.01, confidence=0.7, minlen=2,
               maxlen=4)

fit <- apriori(x, parameter=params)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# remove rhs restrictions again (e)
params <- list(support=0.0075, confidence=0.75,
               minlen=2, maxlen=4)

fit <- apriori(x, parameter=params)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))


