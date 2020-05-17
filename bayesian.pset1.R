## chapter 1, problem 9

#initializations for storage to find median / interval
sims <- 100
i.list <- c(rep(0,sims))
total.num.patients.list <- c(rep(0,sims))
patients.who.waited.list <- c(rep(0,sims))
virtual.line <- c(NA)
time.waited.list <- c(rep(0,sims))
i <- 0

#simulation
for (run in (1:sims)) {
  #initializations
  i <- 0 #time (minutes)
  num.wait <- 0 #the number of patients waiting
  busy <- c(0,0,0) #are or arent the doctors available, 0 if avail, 1 if busy
  free <- c(0,0,0) #when each of the doctors will be free again
  total.num.patients <- 0 #counter for total number of patients who visit
  patients.who.waited <- 0 #counter for number of patients who had to wait
  time.waited <- c(0)
  
  # from 9am to 4pm
  for (i in (1:420)) {
    #if a doctor finishs with a patient at time i, then set corresponding busy to 0 for free
    busy[free==i] <- 0
    #if a doctor finishs with a patient at time i
    if (busy[1]==0 || busy[2]==0 || busy[3]==0) {
      #count the amount of time patient has been waiting ... 
      time.waited[i] <- virtual.line[which(!is.na(virtual.line))[1]]
      #... and then remove from the line
      virtual.line[which(!is.na(virtual.line))[1]] <- NA
    }
    #does a new patient arrive? if yes add one to the waiting line
    if (rexp(1, 0.1)<=1) {
      num.wait <- num.wait+1
      total.num.patients <- total.num.patients+1
      #and add patient to virtual line
      virtual.line[i] <- 1
      #check if patient will have to wait
      if (length(busy[which(busy==0)]) < num.wait){
        patients.who.waited <- patients.who.waited+1
      }
    }
    #check which / how many of the doctors ...
    for (j in (1:3)) {
      #... are not busy, AND check if there is a patient waiting
      if (busy[j]==0 && num.wait > 0){
        #take a patient off the waiting line
        num.wait <- num.wait-1
        #note that doctor is now busy
        busy[j] <- 1
        #set the time when the doctor will no longer be busy to i + pulled value
        free[j] <- i + round(runif(1,15,20),0)
      }
    }
    #if no doctor is free, mark that patients in line had to wait
    for (mark in (which(!is.na(virtual.line)))) {
      virtual.line[mark] <- virtual.line[mark] +1
      }
  }
  # post 4pm
  while (num.wait>0) {
    #if a doctor finishs with a patient at time i, then set corresponding busy to 0
    busy[free==i] <- 0
    #if a doctor finishs with a patient at time i
    if (busy[1]==0 || busy[2]==0 || busy[3]==0) {
      #count the amount of time patient has been waiting ... 
      time.waited[i] <- virtual.line[which(!is.na(virtual.line))[1]]
      #... and then remove from the line
      virtual.line[which(!is.na(virtual.line))[1]] <- NA
    }
    #check which / how many of the doctors ... 
    for (j in (1:3)) {
      #... are not busy, AND check if there is a patient waiting
      if (busy[j]==0 && num.wait > 0){
        #take a patient off the waiting line
        num.wait <- num.wait-1
        #note that that doctor is now busy
        busy[j] <- 1
        #set the time when the doctor will no longer be busy to i + pulled value
        free[j] <- i + round(runif(1,15,20),0)
      }
    }
    #if no doctor is free, mark that patients in line had to wait
    for (mark in (which(!is.na(virtual.line)))) {
      virtual.line[mark] <- virtual.line[mark] +1
    }
    #count elapsed minute
    i <- i + 1
    virtual.line[i] <- NA
  }
  i.list[run] <- i
  total.num.patients.list[run] <- total.num.patients
  patients.who.waited.list[run] <- patients.who.waited
  time.waited.list[run] <- mean(time.waited,na.rm=T)
}

#How many patients came to the office?
median(total.num.patients.list)
total.num.patients.list.ordered <- (sort(total.num.patients.list))
total.num.patients.interval <- total.num.patients.list.ordered[c(25,75)]
total.num.patients.interval

#How many had to wait for a doctor?
median(patients.who.waited.list)
patients.who.waited.list.ordered <- (sort(patients.who.waited.list))
patients.who.waited.list.interval <- patients.who.waited.list.ordered[c(25,75)]
patients.who.waited.list.interval

#What was their average wait?
median(time.waited.list)
time.waited.list.ordered <- (sort(time.waited.list))
time.waited.list.interval <- time.waited.list.ordered[c(25,75)]
time.waited.list.interval

#When did the office close?
#i.list <- i.list-420
median(i.list)
i.list.ordered <- (sort(i.list))
i.list.interval <- i.list.ordered[c(25,75)]
i.list.interval



## chapter 2 problem 4
y <- seq(50,300,1)
dens <- function (x, theta){dnorm (x, 1000*theta, sqrt(1000*theta*(1-theta)))}
dens.mix <- 0.25*dens(y,1/12) + 0.5*dens(y,1/6) + 0.25*dens(y,1/4)
plot (y, dens.mix, ylim=c(0,1.1*max(dens.mix)),
      type="l", xlab="y", ylab="", xaxs="i",
      yaxs="i", yaxt="n", bty="n", cex=2)



##chapter 2 problem 13

# part a
ytilde <- rpois(10000,rgamma(10000,238)/10)
ytilde.ordered <- (sort(ytilde))
interval <- ytilde.ordered[c(250,9750)]
interval

# part b
ytilde <- rpois(10000,rgamma(10000,238)/5.716e12*8e11)
ytilde.ordered <- (sort(ytilde))
interval <- ytilde.ordered[c(250,9750)]
interval

# part c
ytilde <- rpois(10000,rgamma(10000,6919)/10)
ytilde.ordered <- (sort(ytilde))
interval <- ytilde.ordered[c(250,9750)]
interval

# part d
ytilde <- rpois(10000,rgamma(10000,6919)/5.716e12*8e11)
ytilde.ordered <- (sort(ytilde))
interval <- ytilde.ordered[c(250,9750)]
interval
