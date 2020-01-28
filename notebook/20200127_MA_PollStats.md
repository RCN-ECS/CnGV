# Proposed 6 for Poll data

We asked a series of questions and asked workshop attendees to a) assess the importance of questions using likert scale and b) rank top 3 most important questions (hereafter refered to as likert data and rank data, respectively). 

## Likert data 

The most appropriate test for the likert data seems to be the chi-squared test, which tests whether a significant association exists between frequencies observed between 2 categories (e.g., question and rank). However, my goal was to determine exactly which questions were most frequently ranked as "very important", which the chi-squared test does not do. Pairwise comparisons and post-hoc analyses do not work on chi-squared tests nor are they particularly appropriate for this non-normal data.

Instead I decided to perform the hypothesis test using Permutations. Permutations may be appropriate for these data since they dont assume an underlying distribution and *should* allow me to determine which questions were most frequently ranked "very important".

For each question, I calculate the proportion of users that ranked each question as most important. Then I shuffle the data and recalculated 
