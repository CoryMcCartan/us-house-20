# DATA PREPARATION 

library(Hmisc)
library(pollstR)
library(tidyverse)
library(lubridate)

election.day = as.Date("2018-11-06")

# Past Presidential approval
president = read.csv("data/past_pres_approval_raw.csv", sep="\t") %>%
    transmute(approval=Approving / 100,
              date=mdy(End.Date),
              president=President,
              party=Party) %>%
    filter(year(date) >= 1970, year(date) <= 2030) %>%
    arrange(date)
write.csv(president, "data/past_pres_approval.csv", row.names=F)

# Current polling
polls = pollster_charts_polls("2018-national-house-race")[["content"]] %>%
    filter(partisanship == "Nonpartisan") %>%
    transmute(dem = Democrat / 100,
              gop = Republican / 100,
              other = Other / 100,
              undecided = Undecided / 100,
              n_resp = observations,
              n_side = round((dem + gop) * n_resp),
              n_dem = round(dem * n_resp),
              pollster = survey_house,
              firm.id = as.numeric(as.factor(pollster)),
              date = ymd(end_date),
              type = recode(sample_subpopulation, `Registered Voters`="RV",
                            `Likely Votesrs`="LV")) %>%
    as.data.frame
# drop missing data
polls = polls[complete.cases(polls),]
# week IDs
start.day = min(polls$date)
end.day = max(polls$date)
n.weeks = ceiling(as.numeric(election.day - start.day) / 7)
polls$week = floor(as.numeric(election.day - polls$date) / 7)
# output
write.csv(polls, "data/current_polls.csv", row.names=F)

# past polling
past.polls = read.csv("data/raw_past_polling.csv") %>%
    filter(type_simple == "House-G", location == "US")
election.dates = unique(mdy(past.polls$electiondate))
lgt = function(x) log(x / (1-x))
past.polls %<>% transmute(year = year,
                          n_resp = samplesize,
                          dem = cand1_pct/(cand1_pct + cand2_pct),
                          logit = lgt(dem),
                          weeks.until = round(as.numeric(
                              mdy(electiondate) - mdy(polldate)) / 7)) %>%
    group_by(year, weeks.until) %>%
    summarise(lgt = wtd.mean(logit, sqrt(n_resp)),
              sd = sqrt(wtd.var(logit, sqrt(n_resp)))) %>%
    filter(!is.na(sd), sd > 0) %>%
    as.data.frame
# manually add more data
new.pp = data.frame(year=c(1974, 1982, 1994), weeks.until=0,
                    lgt=c(lgt(0.56/0.86), c(0.54/0.9), 0), sd=NA)
election.dates = c(election.dates, ymd("1974-11-01", "1982-11-01", "1994-11-01"))
past.polls = rbind(past.polls, new.pp)

# Merge
economy = read.csv("data/economy.csv", na.strings=".") %>%
    transmute(date = ymd(DATE),
              unemp = as.numeric(UNRATE) / 100,
              gdp = as.numeric(A939RX0Q048SBEA_PC1) / 100,
              infl = as.numeric(CPIAUCSL_PC1) / 100,
              earn = as.numeric(AHETPI_PC1) / 100) %>%
    filter(year(date) > 1970)
congress = read.csv("data/congress_approval.csv")
control = read.csv("data/party_control.csv") %>% 
    mutate(midterm = ifelse((elect.year %% 4) == 2, 1, 0))

fundamentals = data.frame(year=year(election.dates), approval=0, unemp=0, gdp=0, earn=0)
for (i in 1:length(election.dates)) {
    d = election.dates[i]
    appr = (president %>% 
        filter(abs(date - d) < 120 & abs(date - d) > 30) %>%
        summarise(m = mean(approval)))$m
    econ = economy %>%
        filter(abs(date - d) < 120 & abs(date - d) > 30) %>%
        summarise(unemp = mean(unemp),
                  gdp = mean(gdp),
                  earn = mean(earn))
    
    fundamentals[i,2] = log(appr / (1 - appr))
    fundamentals[i,3] = econ$unemp
    fundamentals[i,4] = econ$gdp
    fundamentals[i,5] = econ$earn
}

code = function(x) 2*x - 1

combined = past.polls %>%
    left_join(control, by=c("year"="elect.year")) %>%
    left_join(fundamentals, by=c("year")) %>%
    transmute(year, weeks_until=weeks.until, logit_intent=lgt, sd_intent=sd,
           seats=dem.seats, before=dem.before, pres=code(dem.pres), 
           house=code(dem.house), midterm=midterm, appr=approval, unemp, gdp, earn)

# manually add 1974 and fix 2014
combined$weeks_until[combined$year==2014][1] = 0
combined[combined$year==1974, c("seats", "before", "pres", "house", "midterm")] =
    c(291, 242, -1, 1, 1)

write.csv(combined, "data/combined.csv", row.names=F)
