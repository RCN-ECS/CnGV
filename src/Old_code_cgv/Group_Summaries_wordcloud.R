library(wordcloud)
library(SnowballC)
library(tm)
library(RColorBrewer)

GSum<-readLines("~/Desktop/Group_Summaries_All.txt")

docs <- Corpus(VectorSource(GSum))
inspect(docs)

toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
docs <- tm_map(docs, toSpace, "-")
docs <- tm_map(docs, toSpace, "'")
docs <- tm_map(docs, toSpace, "\\|")

toSingle <- content_transformer(function (x , pattern ) gsub('.{1}$', '', x)) #Function to remove last letter of word
noDash <- content_transformer(function (x , pattern ) gsub('-', ' ', x)) #Function to hyphens
noAsterisk <- content_transformer(function (x , pattern ) gsub('*', ' ', x)) #Function to remove asterisks

docs <- tm_map(docs, toSingle, "populations") #Change "populations" to "population"
docs <- tm_map(docs, toSingle, "genomics") #Change "genomics" to "genomic"
docs <- tm_map(docs, toSingle, "traits") #Change "traits" to "trait"
docs <- tm_map(docs, noDash) #Remove dash
docs <- tm_map(docs, noAsterisk) #Remove asterisks
docs <- tm_map(docs, content_transformer(tolower)) #Change all letters to lower case 
docs <- tm_map(docs, removeNumbers) #Remove all numbers
docs <- tm_map(docs, removeWords, stopwords("english")) #Remove common english stop words
docs <- tm_map(docs, removeWords, c("the", "and", "that", "for", "are", "with","may", 
                                    "can", "will","but","not","more","could","this",
                                    "how", "which","would","have", "when","there",
                                    "then", "from", "what", "should", "they")) #Remove commonly encountered non-interesting words
docs <- tm_map(docs, removePunctuation) #Remove punctuation
docs <- tm_map(docs, stripWhitespace) #Strip away white space
#docs <- tm_map(docs, stemDocument) #Text stemming


dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 10)


set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=500, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Most frequent words",
        ylab = "Word frequencies")

