setwd("~/Documents/data/CapstoneProject/final/en_US")

library(tm)
library(stringr)
library(dplyr)
library(ggplot2)

# set options
options(stringsAsFactors = FALSE)

# get list of files
flist <- list.files()

# define a function to sample from each files
file_sample <- function(f) {
  frac = 1 # fraction of lines to use in the model
  lines <- readLines(f, skipNul = TRUE)
  smpl <- sample(lines, round(frac*length(lines)), replace = FALSE)
  #smpl <- lines   # we use all the lines here
  return(c(f,smpl))
}

# apply the function to all the files 
samples_list <- lapply(flist, file_sample)

# create a "Volatile Corpus" by putting together samples from all 3 files
docs <- VCorpus(VectorSource(samples_list))

# assign ID to each part of the corpus
for (i in 1:length(docs)) {
  attr(docs[[i]], "ID") <- samples_list[[i]][1]
}

# define necessary functions for pre-processing
removeForeign <- function(x) iconv(x,"latin1","ASCII",sub = "")
removeURL <- function(x) gsub("http:[[:alnum:]]*", "", x)
removeHashtag <- function(x) gsub("#\\S+", "", x)
removeUsername <- function(x) gsub("@\\S+", "", x)

PreProcess <- function(doc){
  #doc <- tm_map(doc, removeWords, cursewords)
  doc <- tm_map(doc, removeWords, stopwords("english"))
  doc <- tm_map(doc, removePunctuation)
  doc <- tm_map(doc, removeNumbers)
  doc <- tm_map(doc, stripWhitespace)
  doc <- tm_map(doc, content_transformer(tolower))
  doc <- tm_map(doc, content_transformer(removeForeign))
  doc <- tm_map(doc, content_transformer(removeURL))
  doc <- tm_map(doc, content_transformer(removeHashtag))
  doc <- tm_map(doc, content_transformer(removeUsername))
  #doc <- tm_map(doc, PlainTextDocument)
  return(doc)
}

# apply the pre-processing function on the corpus
docs <- PreProcess(docs)

# define required functions for tokenization
BigramTokenizer <- function(x) unlist(lapply(ngrams(words(x), 2), paste, collapse = " "), use.names = FALSE)
TrigramTokenizer <- function(x) unlist(lapply(ngrams(words(x), 3), paste, collapse = " "), use.names = FALSE)

# tokenize the text
Unigram_tdm <- TermDocumentMatrix(docs)
Bigram_tdm <- TermDocumentMatrix(docs, control = list(tokenize = BigramTokenizer))
Trigram_tdm <- TermDocumentMatrix(docs, control = list(tokenize = TrigramTokenizer))


# number of most frequent n-grams to show
n <- 20

# find most frequent n-grams
UnigramFreq <- findFreqTerms(Unigram_tdm, n)
BigramFreq <- findFreqTerms(Bigram_tdm, n)
TrigramFreq <- findFreqTerms(Trigram_tdm, n)

# generate sorted data frames
UnigramFrequencies <- rowSums(as.matrix(Unigram_tdm[UnigramFreq,]))
UnigramFrequencies <- data.frame(ngram=names(UnigramFrequencies),Frequency=UnigramFrequencies)
UnigramFrequencies <- UnigramFrequencies[order(-UnigramFrequencies$Frequency),]

BigramFrequencies <- rowSums(as.matrix(Bigram_tdm[BigramFreq,]))
BigramFrequencies <- data.frame(ngram=names(BigramFrequencies),Frequency=BigramFrequencies)
BigramFrequencies <- BigramFrequencies[order(-BigramFrequencies$Frequency),]

TrigramFrequencies <- rowSums(as.matrix(Trigram_tdm[TrigramFreq,]))
TrigramFrequencies <- data.frame(ngram=names(TrigramFrequencies),Frequency=TrigramFrequencies)
TrigramFrequencies <- TrigramFrequencies[order(-TrigramFrequencies$Frequency),]


# set bigram and trigram discounts
gamma2 <- 0.5  # bigram discount
gamma3 <- 0.5  # trigram discount

# bigram prefix of word we want to predict
bigPre <- 'on_my'

# estimate probabilities of words completing the observed 3-grams

getObsTrigs <- function(bigPre, trigs) {
  trigs.winA <- data.frame(ngram=vector(mode = 'character', length = 0),
                           Frequency=vector(mode = 'integer', length = 0))
  regex <- sprintf("%s%s%s", "^", bigPre, "_")
  trigram_indices <- grep(regex, trigs$ngram)
  if(length(trigram_indices) > 0) {
    trigs.winA <- trigs[trigram_indices, ]
  }
  
  return(trigs.winA)
}

getObsTriProbs <- function(obsTrigs, bigrs, bigPre, triDisc=0.5) {
  if(nrow(obsTrigs) < 1) return(NULL)
  obsCount <- filter(bigrs, ngram==bigPre)$Frequency[1]
  obsTrigProbs <- mutate(obsTrigs, Frequency=((Frequency - triDisc) / obsCount))
  colnames(obsTrigProbs) <- c("ngram", "probability")
  
  return(obsTrigProbs)
}

obs_trigs <- getObsTrigs(bigPre, TrigramFrequencies)  # get trigrams and counts

# convert counts to probabilities
qbo_obs_trigrams <- getObsTriProbs(obs_trigs, BigramFrequencies, bigPre, gamma3)
qbo_obs_trigrams

getUnobsTrigTails <- function(obsTrigs, unigs) {
  obs_trig_tails <- str_split_fixed(obsTrigs, "_", 3)[, 3]
  unobs_trig_tails <- unigs[!(unigs$ngram %in% obs_trig_tails), ]$ngram
  return(unobs_trig_tails)
}

unobs_trig_tails <- getUnobsTrigTails(obs_trigs$ngram, UnigramFrequencies)
unobs_trig_tails


getAlphaBigram <- function(unigram, bigrams, bigDisc=0.5) {
  # get all bigrams that start with unigram
  regex <- sprintf("%s%s%s", "^", unigram$ngram[1], "_")
  bigsThatStartWithUnig <- bigrams[grep(regex, bigrams$ngram),]
  if(nrow(bigsThatStartWithUnig) < 1) return(0)
  alphaBi <- 1 - (sum(bigsThatStartWithUnig$Frequency - bigDisc) / unigram$Frequency)
  
  return(alphaBi)
}

unig <- str_split(bigPre, "_")[[1]][2]
unig <- UnigramFrequencies[UnigramFrequencies$ngram == unig,]
alpha_big <- getAlphaBigram(unig, BigramFrequencies, gamma2)
alpha_big


getBoBigrams <- function(bigPre, unobsTrigTails) {
  w_i_minus1 <- str_split(bigPre, "_")[[1]][2]
  boBigrams <- paste(w_i_minus1, unobsTrigTails, sep = "_")
  return(boBigrams)
}


getObsBoBigrams <- function(bigPre, unobsTrigTails, bigrs) {
  boBigrams <- getBoBigrams(bigPre, unobsTrigTails)
  obs_bo_bigrams <- bigrs[bigrs$ngram %in% boBigrams, ]
  return(obs_bo_bigrams)
}

getUnobsBoBigrams <- function(bigPre, unobsTrigTails, obsBoBigram) {
  boBigrams <- getBoBigrams(bigPre, unobsTrigTails)
  unobs_bigs <- boBigrams[!(boBigrams %in% obsBoBigram$ngram)]
  return(unobs_bigs)
}

getObsBigProbs <- function(obsBoBigrams, unigs, bigDisc=0.5) {
  first_words <- str_split_fixed(obsBoBigrams$ngram, "_", 2)[, 1]
  first_word_freqs <- unigs[unigs$ngram %in% first_words, ]
  obsBigProbs <- (obsBoBigrams$Frequency - bigDisc) / first_word_freqs$Frequency
  obsBigProbs <- data.frame(ngram=obsBoBigrams$ngram, probability=obsBigProbs)
  
  return(obsBigProbs)
}

getQboUnobsBigrams <- function(unobsBoBigrams, unigs, alphaBig) {
  # get the unobserved bigram tails
  qboUnobsBigs <- str_split_fixed(unobsBoBigrams, "_", 2)[, 2]
  w_in_Aw_iminus1 <- unigs[!(unigs$ngram %in% qboUnobsBigs), ]
  # convert to data.frame with counts
  qboUnobsBigs <- unigs[unigs$ngram %in% qboUnobsBigs, ]
  denom <- sum(qboUnobsBigs$Frequency)
  # converts counts to probabilities
  qboUnobsBigs <- data.frame(ngram=unobsBoBigrams,
                             probability=(alphaBig * qboUnobsBigs$Frequency / denom))
  
  return(qboUnobsBigs)
}

bo_bigrams <- getBoBigrams(bigPre, unobs_trig_tails)  # get backed off bigrams

obs_bo_bigrams <- getObsBoBigrams(bigPre, unobs_trig_tails, BigramFrequencies)
unobs_bo_bigrams <- getUnobsBoBigrams(bigPre, unobs_trig_tails, obs_bo_bigrams)

qbo_obs_bigrams <- getObsBigProbs(obs_bo_bigrams, UnigramFrequencies, gamma2) #ngram     probs

unig <- str_split(bigPre, "_")[[1]][2]
unig <- UnigramFrequencies[UnigramFrequencies$ngram == unig,]

qbo_unobs_bigrams <- getQboUnobsBigrams(unobs_bo_bigrams, UnigramFrequencies, alpha_big)
qbo_bigrams <- rbind(qbo_obs_bigrams, qbo_unobs_bigrams)
qbo_bigrams

getAlphaTrigram <- function(obsTrigs, bigram, triDisc=0.5) {
  if(nrow(obsTrigs) < 1) return(1)
  alphaTri <- 1 - sum((obsTrigs$Frequency - triDisc) / bigram$Frequency[1])
  
  return(alphaTri)
}

bigram <- BigramFrequencies[BigramFrequencies$ngram %in% bigPre, ]
alpha_trig <- getAlphaTrigram(obs_trigs, bigram, gamma3)
alpha_trig


getUnobsTriProbs <- function(bigPre, qboObsBigrams,
                             qboUnobsBigrams, alphaTrig) {
  qboBigrams <- rbind(qboObsBigrams, qboUnobsBigrams)
  qboBigrams <- qboBigrams[order(-qboBigrams$probability), ]
  sumQboBigs <- sum(qboBigrams$probability)
  first_bigPre_word <- str_split(bigPre, "_")[[1]][1]
  unobsTrigNgrams <- paste(first_bigPre_word, qboBigrams$ngram, sep="_")
  unobsTrigProbs <- alphaTrig * qboBigrams$probability / sumQboBigs
  unobsTrigDf <- data.frame(ngram=unobsTrigNgrams, probability=unobsTrigProbs)
  
  return(unobsTrigDf)
}

qbo_unobs_trigrams <- getUnobsTriProbs(bigPre, qbo_obs_bigrams,
                                       qbo_unobs_bigrams, alpha_trig)
qbo_unobs_trigrams

getPredictionMsg <- function(qbo_trigs) {
  # pull off tail word of highest prob trigram
  prediction <- str_split(qbo_trigs$ngram[1], "_")[[1]][3]
  result <- sprintf("%s%s%s%.4f", "highest prob prediction is >>> ", prediction,
                    " <<< which has probability = ", qbo_trigs$probability[1])
  return(result)
}

qbo_trigrams <- rbind(qbo_obs_trigrams, qbo_unobs_trigrams)
qbo_trigrams <- qbo_trigrams[order(-qbo_trigrams$probability), ]  # sort by desc prob
out_msg <- getPredictionMsg(qbo_trigrams)
out_msg

