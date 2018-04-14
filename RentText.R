### TEXT ANALYSIS ######
# here we will do an exploratory analysis of some real nyc text data


#install.packages("tidytext")
#install.packages("dplyr")
#install.packages("stringr")

library(tidytext)
library(dplyr)
library(stringr)

# set the working directory to wherever you've downloaded all
# the files
setwd("C:/Users/owenw/Documents/jla2167_ddseminar/jla2167_ddseminar/dat/txt/")
text_df <- data.frame()

for(file in dir()){
  dat <- substr(file,1,10)  
  temp <- readChar(file,file.info(file)$size)
  rm_string <- function(pattern, x) gsub(x = x, pattern = pattern , replacement = " " )
  rmv <- c("[\r\n\t]","[0-9]+",
           "ROY ALLEN & ASSOCIATES","ROY ALLEN &   ASSOCIATES, INC.",
           "INC","521 Fifth Avenue","17th Floor","New York, New York")
  temp1 <- Reduce(rm_string, rmv, init = temp, right = T) 
  text <- unlist(strsplit(temp1,':',fixed=TRUE))
  people <- sub(patt='(.+)(([ ,.]+\\w+){2})[ ]?$', repl='\\2', text)
  people <- word(people,-1)[-1]
  text <- rm_string("\\w*$",trimws(rm_string("\\w*$",text)))
  rmv <- c("[[:punct:]]","MR.","MS.")
  text <- trimws(Reduce(rm_string, rmv, init = text, right = T))
  text <- gsub("   "," ",text) 
  text <- gsub("  "," ",text) 
  text <- text[toupper(people) == people]
  people <- people[toupper(people) == people]
  people <- people[nchar(text) > 50]
  text <- text[nchar(text) > 50]
  
  text <- text[-c(1:5,(length(people)-5):length(people))]
  people <- people[-c(1:5,(length(people)-5):length(people))]
  
  text_df <-  rbind(df,data.frame(date=rep(dat,length(people)),
                             people=people,
                             test=text)
  )
}

# dataset quite large so may be good to test stuff on just small
# sample
#text_df <- text_df[1:10,]
# this extracts the individual words from all the text
text_df$test = as.character(text_df$test)
text_df_words <- text_df %>%
  unnest_tokens(word, test )

text_df_words

# some common stop words which may not be useful
data(stop_words)
stop_words

# this removes the stop words from the data
text_df_main <- text_df_words %>%
  anti_join(stop_words)

# counts the most common words in the data
text_df_words %>%
  count(word, sort = TRUE)

text_df_main_count<- text_df_main %>%
  count(word,sort = TRUE)

# print the 100 most common words
print(text_df_main_count,n=100)


library(wordcloud)

# might want to add custom stop words here 
# like the names, director, board, etc
custom_stop_words <- bind_rows(data_frame(word = c("director"), 
                                          lexicon = c("custom")), 
                               stop_words)

text_df_words %>%
  anti_join(stop_words) %>% # use custom stop words here
  count(word) %>%
  with(wordcloud(word, n, max.words = 75))

# Try some sentiment analysis ####
# try to evaluate which words are positive or negative
sentiments
get_sentiments("bing")

word_sentiments <-text_df_words %>%
  anti_join(custom_stop_words) %>%
  inner_join(get_sentiments("bing")) %>%
  count(word, sentiment,sort = TRUE) %>%
  ungroup()

library(ggplot2)

word_sentiments %>%
  group_by(sentiment) %>%
  top_n(10) %>%
  ungroup() %>%
  mutate(word = reorder(word, n)) %>%
  ggplot(aes(word, n, fill = sentiment)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~sentiment, scales = "free_y") +
  labs(y = "Contribution to sentiment",
       x = NULL) +
  coord_flip()

library(reshape2)

text_df_words %>%
  inner_join(get_sentiments("bing")) %>%
  count(word, sentiment, sort = TRUE) %>%
  acast(word ~ sentiment, value.var = "n", fill = 0) %>%
  comparison.cloud(colors = c("gray20", "gray80"),
                   max.words = 100)

# try to do sentiment analysis at the sentence level
text_sentences <- text_df %>% 
  unnest_tokens(sentence, test, token = "sentences")

text_sentences[1,]
text_sentences$sentence[1]

# this is clearly more than one sentence! so the data when cleaned must
# have removed ., etc, so would need to clean it differently to fix this


library(dplyr)
library(janeaustenr)
library(tidytext)

book_words <- austen_books() %>%
  unnest_tokens(word, text) %>%
  count(book, word, sort = TRUE) %>%
  ungroup()

total_words <- book_words %>% 
  group_by(book) %>% 
  summarize(total = sum(n))

book_words <- left_join(book_words, total_words)

book_words

# Zipf's law states that the frequency that a word appears is inversely proportional to its rank. 

library(ggplot2)

ggplot(book_words, aes(n/total, fill = book)) +
  geom_histogram(show.legend = FALSE) +
  xlim(NA, 0.0009) +
  facet_wrap(~book, ncol = 2, scales = "free_y")

freq_by_rank <- book_words %>% 
  group_by(book) %>% 
  mutate(rank = row_number(), 
         `term frequency` = n/total)

freq_by_rank

freq_by_rank %>% 
  ggplot(aes(rank, `term frequency`, color = book)) + 
  geom_line(size = 1.1, alpha = 0.8, show.legend = FALSE) + 
  scale_x_log10() +
  scale_y_log10()

rank_subset <- freq_by_rank %>% 
  filter(rank < 500,
         rank > 10)

lm(log10(`term frequency`) ~ log10(rank), data = rank_subset)

# freq prop to 1/rank would be slope of -1 so pretty close here

freq_by_rank %>% 
  ggplot(aes(rank, `term frequency`, color = book)) + 
  geom_abline(intercept = -0.62, slope = -1.1, color = "gray50", linetype = 2) +
  geom_line(size = 1.1, alpha = 0.8, show.legend = FALSE) + 
  scale_x_log10() +
  scale_y_log10()


library(topicmodels)

data("AssociatedPress")
AssociatedPress

# set a seed so that the output of the model is predictable
ap_lda <- LDA(AssociatedPress, k = 2, control = list(seed = 1234))
ap_lda

library(tidytext)

ap_topics <- tidy(ap_lda, matrix = "beta")
ap_topics

library(ggplot2)
library(dplyr)

ap_top_terms <- ap_topics %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

ap_top_terms %>%
  mutate(term = reorder(term, beta)) %>%
  ggplot(aes(term, beta, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  coord_flip()

library(tidyr)

beta_spread <- ap_topics %>%
  mutate(topic = paste0("topic", topic)) %>%
  spread(topic, beta) %>%
  filter(topic1 > .001 | topic2 > .001) %>%
  mutate(log_ratio = log2(topic2 / topic1))

beta_spread


