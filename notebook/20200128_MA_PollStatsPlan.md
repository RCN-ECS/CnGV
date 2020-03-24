# Possible Stats Plan for Poll Results

I am seeking the most appropriate statistical method for the consortium data (with the full acknowledgement that stats do not need to be fancy or heavy handed in this situation).

We asked two complementary families of questions (other than demographic info). We asked workshop attendees to rank each question that emerged from the workshop according to a 3 point likert scale (hereafter called Likert data), and we asked attendees to rank the top 3 questions with the greatest ability to push the fields of marine science, evolutionary biology, and conservation forward (hereafter refered to as rank data).

## Likert Data

For the likert data, Chi-squared is probably the most appropriate but doesnt allow or post-hoc or pairwise comparisons which I what I want. 
I tried a different approach using permutations. For each question, I calculated the proportion of users that ranked the question as "Very Important". I then reshuffled all responses among all categories. After shuffling, I re-calculated the proportion ranked as "Very Important" and repeated this 999 times. 

I then use that null distribution to determine for each question, whether the proportion of individuals that ranked the question as very important was significant. 

Importantly, this approach measures significance by comparing each question within itself, and does NOT compare across questions. I'd prefer to compare across questions, but I'm not sure how to do that at this moment. 

When divided, the Marine Science and Evolutionary biology categories didn't have any that were significant at alpha = 0.05. I plotted the effect size (percent ranked) against the pvalue. 

Perhaps a way to choose the "top" questions are to choose the questions with the highest rank percent with lowest p-value.

Below are the three plots showing the results of the permutation. 
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/MarineScience_Permuteplot.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/Evolution_permutePlot.png)
![image](https://github.com/RCN-ECS/CnGV/blob/master/results/MarineEvo_PermutePlot.png)

## Rank Data

I haven't worked with this yet. If this permutation approach is okay, I'll do something similar. 

For each question I'll estimate the proportion of users that ranked it. Then reshuffle to find out which questions were most frequently ranked as most important.

