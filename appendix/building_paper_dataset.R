##############################################################
## Script: Scrapping Papers for Distributive Politics Paper ##
## Authors: Alptekin, Freire, Mignozzetti, and Roman        ##
## Last Modified: Mar 12, 2021                              ##
##############################################################

# Note: this script is out of order, and contains
#       the code for scraping the article's data

# Required packages
pkgs <- c("tidyverse", "rvest", "RSelenium")

# Install the packages if necessary
installIfNot <- function(x) {
  if (x %in% rownames(installed.packages()) == FALSE)
    install.packages(x, dependencies = T,
                     repos = "http://cran.us.r-project.org")
}
lapply(pkgs, installIfNot)

# Load packages
lapply(pkgs, require, character.only = T)
rm(pkgs, installIfNot)

# Setting Up Selenium

# Alternative 1: Setting up Selenium (head)
rsD <- rsDriver(port = 1114L, browser = c("firefox"))
remDr <- rsD$client
remDr$open()

# Google Scholar

# site: https://scholar.google.com
remDr$navigate("https://scholar.google.com/
               scholar?cites=13117579863846712459&as_sdt=2005&sciodt=0,5")

articles_weingast <- tibble(
  value = NA,
  term = NA,
  page = NA
)

k <- 0

for (j in 1:213) { # we had to manually choose the number of pages here
  
  Sys.sleep(rpois(1, 5))
  
  # Getting articles basic information
  k <- k + 1
  
  webElem <- remDr$findElement("css", "body")
  
  
  title <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(
      xpath = '//*[contains(concat( " ", @class, " " ),
      concat( " ", "gs_rt", " " ))]'
    ) %>%
    html_text() %>%
    enframe(name = NULL) %>%
    rename("title" = "value") %>%
    mutate(page = k)
  
  articles_partial <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(xpath = '//*[contains(concat( " ", @class, " " ), concat( " ", "gs_a", " " ))]') %>%
    html_text() %>%
    enframe(name = NULL) %>%
    bind_rows(
      tibble(
        value = "delete",
        term = NA,
        page = NA
      ),
      .
    ) %>%
    bind_cols(., title)
  
  
  # Binding articles
  articles_weingast <- bind_rows(articles_weingast, articles_partial)
  
  # Changing Pages
  next_button <- remDr$findElement(using = "xpath", "/html/body/div/div[11]/div[2]/div[2]/div[3]/
                                   div[2]/center/table/tbody/tr/td[12]/a/b")
  next_button$clickElement()
  # Deleting cookies
  remDr$deleteAllCookies()
}

write_csv(articles_weingast, "scholar_weingast_raw.csv")

articles_weingast <- articles_weingast %>%
  select(-term, -page) %>%
  filter(value != "delete") %>%
  slice(2:nrow(.)) %>%
  filter(!grepl("books.google.com", value)) %>%
  filter(!grepl("BOOK", title)) %>%
  separate(., value, into = c("author", "value"),
           sep = " -", remove = T, extra = "merge", fill = "right") %>%
  separate(., value, into = c("journal", "year"),
           sep = ",", remove = T, extra = "merge", fill = "right") %>%
  mutate(
    year_2 = ifelse(is.na(year), journal, year),
    journal = ifelse(is.na(year), NA, journal),
    year = year_2
  ) %>%
  select(-year_2) %>%
  separate(., year, into = c("year", "site"),
           sep = "-", remove = T, extra = "merge", fill = "right") %>%
  mutate(
    year = gsub("[^0-9 ]", "", value),
    year = gsub("^[0-9]{5,}", "", year),
    year = gsub("^ {1,}", "", year)
  ) %>%
  separate(., year, into = c("year", "junk"),
           sep = " ", remove = T, extra = "merge", fill = "right") %>%
  separate(., value, into = c("journal", "value"),
           sep = "-", remove = T, extra = "merge", fill = "right") %>%
  mutate(
    journal_untidy = year,
    year = gsub("[^0-9 ]", "", year)
  )

articles <- articles %>%
  na.omit() %>%
  distinct(., value, .keep_all = T)

write_csv(articles, "google_scholar_clean.csv")

# Scopus

# url: https://www.scopus.com/home.uri

# Scraping Scopus requires a bit more manual labor.
# You can login on scopus through your university/institution
# and download the metadata of the article(s) you
# want directly from there.
# All we need to do after that is scrape the information of
# every link from the .csv file downloaded previously

scopus <- read_csv("scopus.csv")

scopus <- scopus %>%
  mutate(article = map_chr(Link, ~ {
    remDr$navigate(.x)
    read_html(remDr$getPageSource()[[1]]) %>%
      html_nodes(xpath = '//*[(@id = "abstractSection")]//p') %>%
      html_text() %>%
      paste(., collapse = "\r\n")
  })) %>%
  mutate(
    article = gsub(
      '\r\n\nUse this section.*Topics\n\n\n"',
      "",
      article
    ),
    article = gsub(
      "Topics are unique.*onwards.",
      "",
      article
    ),
    article = gsub(
      "Use this section.*documents.",
      "",
      article
    ),
    article = gsub(
      "Learn more about these Topics",
      "",
      article
    ),
    article = gsub(
      ". 20.*, Springer Science\\+Business Media, LLC, part of Springer Nature.",
      "",
      article
    ),
    article = gsub("\r", "", article),
    article = gsub("\n", "", article),
    article = gsub(" {2,}", " ", article),
    article = gsub(" {3,}", "", article)
  )

write_csv(scopus, "scopus_clean.csv")

#### Microsoft Academic ####

# url: https://academic.microsoft.com/home
articles <- list()
k <- 0
remDr$navigate("https://academic.microsoft.com/paper/
               2076316673/citedby/search?q=The%20Political%
               20Economy%20of%20Benefits%20and%20Costs%3A%20A%
               20Neoclassical%20Approach%20to
               %20Distributive%20Politics&qe=RId%253D2076316673&f=&orderBy=0")


# 1) Getting hyperlinks from articles
for (j in 1:100) {
  k <- k + 1
  print(k)
  
  # Navigating Website
  
  Sys.sleep(rpois(2, 5))
  
  # Getting articles' links
  articles[[k]] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(xpath = "//a") %>%
    html_attr("href") %>%
    enframe(name = NULL) %>%
    filter(
      grepl("paper/", value),
      !grepl("citedby", value)
    ) %>%
    mutate(value = paste0("https://academic.microsoft.com/", value))
  
  # Changing Page
  next_page_others <- remDr$findElement(using = "xpath", "/html/body/div/div/div/router-view/router-view/                                        ma-edp-serp/div/div[2]/div/
  compose/div/div[2]/ma-pager/div/i[2]")
  next_page_others$clickElement()
}


articles <- articles %>%
  reduce(bind_rows) %>%
  distinct(value, .keep_all = T)

# 2) Navigating through articles and scraping them

articles_links <- articles %>%
  mutate(
    abstract = NA,
    title = NA,
    year = NA,
    journal = NA,
    authors = NA,
    tags = NA
  )

for (i in 1:nrow(articles_links)) {
  remDr$navigate(articles_links$value[i])
  Sys.sleep(rpois(1, 4))
  
  articles_links$abstract[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(., xpath = "//html/body/div/div/div/router-view/compose[1]/
               div/div/ma-entity-detail-info/compose/div/div/div[1]/p") %>%
    html_text() %>%
    paste(., collapse = " ")
  
  articles_links$title[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(.,
               xpath = '//*[contains(concat( " ", @class, " " ),
               concat( " ", "name", " " ))]') %>%
    html_text() %>%
    paste(., collapse = " ")
  
  articles_links$year[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(.,
               xpath = '//*[contains(concat( " ", @class, " " ),
               concat( " ", "name-section", " " ))]
               //*[contains(concat( " ", @class, " " ), concat( " ", "year", " " ))]') %>%
    html_text() %>%
    paste(., collapse = " ")
  
  articles_links$journal[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(., xpath = '//*[contains(concat( " ", @class, " " )
               , concat( " ", "pub-name", " " ))]') %>%
    html_text() %>%
    paste(., collapse = " ")
  
  articles_links$authors[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(., xpath = "/html/body/div/div/div/router-view/compose[1]/div/div/
               ma-entity-detail-info/compose/div/div/div[1]/
               ma-author-string-collection") %>%
    html_text() %>%
    paste(., collapse = " ")
  
  articles_links$tags[i] <- read_html(remDr$getPageSource()[[1]]) %>%
    html_nodes(., xpath = "/html/body/div/div/div/router-view/compose[1]/
               div/div/ma-entity-detail-info/compose/
               div/div/div[1]/ma-tag-cloud/div") %>%
    html_text() %>%
    paste(., collapse = " ")
}

write_csv(articles_links, "microsoft_academic_clean.csv")