library(tidyverse)
library(lubridate)
library(fredr)


# PRESIDENTIAL APPROVAL
president = read_csv("data/past_pres_approval.csv")
curr_pres = read_csv("../house-data/trump_polls.csv") %>%
    transmute(approval = approve/100,
              date = mdy(enddate),
              president = president,
              party = "Republican")
president =  bind_rows(president, curr_pres)

# GENERIC BALLOT POLLING
past_polls = read_csv("../house-data/raw_past_polling.csv") %>%
    filter(type_simple == "House-G", location == "US")
election_dates = unique(mdy(past_polls$electiondate))
past_polls = past_polls %>%
        transmute(year = year,
                  n_resp = samplesize,
                  dem = cand1_pct/(cand1_pct + cand2_pct),
                  logit = qlogis(dem),
                  weeks_until = round(as.numeric(mdy(electiondate) -
                                                     mdy(polldate)) / 7)) %>%
    group_by(year, weeks_until) %>%
    summarise(lgt = weighted.mean(logit, sqrt(n_resp)),
              sd = sqrt(weighted.mean((logit - lgt)^2, sqrt(n_resp))) ) %>%
    filter(!is.na(sd), sd > 0) %>% ungroup
# manually add more data
new_pp = tibble(year=c(1974, 1976, 1978, 1982, 1994), weeks_until=0,
                    lgt=c(qlogis(0.56/0.86), NA, NA, qlogis(0.54/0.9), 0), sd=0.13)
#                                      imputed with max sd previously ^^^^
election_dates = c(election_dates, ymd("1974-11-01", "1976-11-02", "1978-11-07",
                                       "1982-11-01", "1994-11-01"))
past_polls = bind_rows(past_polls, new_pp)


# ECONOMY
fredr_set_key(Sys.getenv("FRED_KEY"))
econ_params = list(
    series = c("UNRATE", "A939RX0Q048SBEA", "CPIAUCSL", "AHETPI"),
    units = c("lin", "pc1", "pc1", "pc1")
)
min_econ_date = ymd("1970-01-01")
economy = pmap_dfr(econ_params, ~ fredr(.x, observation_start=min_econ_date,
                              frequency="q", units=.y)) %>%
    pivot_wider(names_from=series_id, values_from=value) %>%
    transmute(date = date,
              unemp = as.numeric(UNRATE) / 100,
              gdp = as.numeric(A939RX0Q048SBEA) / 100,
              infl = as.numeric(CPIAUCSL) / 100,
              earn = as.numeric(AHETPI) / 100)
write_csv(economy, "data/economy.csv")

congress = read_csv("data/congress_approval.csv")
control = read_csv("data/party_control.csv") %>%
    mutate(midterm = ifelse((elect_year %% 4) == 2, 1, 0))

# COMBINE
election_dates = ymd(str_c(seq(1974, 2018, 2), "-11-01"))
fundamentals = tibble(year=year(election_dates), approval=0, unemp=0, gdp=0, earn=0)
for (i in 1:length(election_dates)) {
    d = election_dates[i]
    appr = (president %>%
        filter(abs(date - d) < 120 & abs(date - d) > 60) %>%
        summarise(m = mean(approval)))$m
    econ = economy %>%
        filter(abs(date - d) < 120 & abs(date - d) > 60) %>%
        summarise(unemp = mean(unemp),
                  gdp = mean(gdp),
                  earn = mean(earn))

    fundamentals[i,2] = log(appr / (1 - appr))
    fundamentals[i,3] = econ$unemp
    fundamentals[i,4] = econ$gdp
    fundamentals[i,5] = econ$earn
}

code = function(x) 2*x - 1

combined = rename(control, year=elect_year) %>%
    #left_join(control, by=c("year"="elect_year")) %>%
    left_join(fundamentals, by=c("year")) %>%
    transmute(year=year, #logit_intent=lgt, sd_intent=sd,
              seats=dem_seats, before=dem_before,
              retire_lgt = qlogis(dem_retire / (gop_retire + dem_retire)),
              logit_pop = log(dem_pop / gop_pop),
              pres=code(dem_pres), house=code(dem_house), midterm=midterm,
              appr=approval, unemp, gdp, earn) %>%
    arrange(year) %>%
    unique

write_csv(combined, "data/combined.csv")
