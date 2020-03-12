library(dplyr)
library(forcats)
library(readr)
library(lubridate)

load_polls = function(start_date=ymd("2019-11-12"),
                      election_day=ymd("2020-11-03"), write=T) {

    # 2020
    polls_url = "https://projects.fivethirtyeight.com/generic-ballot-data/generic_polllist.csv"
    polls_raw = suppressMessages(read_csv(polls_url))

    # week IDs
    start_date = ymd("2019-11-12") # one week after election day
    n_days = as.numeric(election_day - start_date)
    n_weeks = ceiling(as.numeric(election_day - start_date) / 7)

    polls = polls_raw %>%
        mutate(date = mdy(startdate) + (mdy(enddate) - mdy(startdate))/2,
               week = n_weeks - ceiling(as.numeric(election_day - date) / 7) + 1,
               n_weeks = n_weeks,
               n_dem = round(samplesize * dem/100),
               twoparty = dem + rep,
               sample = pmin(round(samplesize * twoparty/100), 3000),
               dem = dem / twoparty,
               gop = rep / twoparty,
               logit_inflate = 1/dem + 1/(1 - dem),
               var_poll = logit_inflate^2 * dem*(1 - dem)/sample * (1+coalesce(tracking, F)),
                        # approx adjust to account for increased tracking error
               firm = as.character(fct_lump(pollster, 30)),
               #firm = abbreviate(pollster, 12, named=F),
               type_rv = population=="rv" | population == "v",
               type_lv = population=="lv",
               type_a = population=="a") %>%
        select(date, week, n_weeks, firm, sample, type_rv, type_lv, type_a,
               dem, n_dem, var_poll) %>%
        filter(date >= start_date, !is.na(var_poll))

    if (write) write_csv(polls, "data/current_polls.csv")

    polls
}
