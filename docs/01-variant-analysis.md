🧬 Track 1 - Variant Prediction
The Goal: Predict the specific genetic variant(s) that may be driving the child's condition using their clinical symptoms and genomic data.
You'll be given the child's genomic data (FASTQ + VCF) and a description of their symptoms, and asked to predict which specific variant(s) may be causing the condition. Secondary or incidental findings are welcome too.

Your ranked predictions are automatically checked against the clinically confirmed answer and your score updates live on the public leaderboard.

What you submit: A ranked list of variant predictions, a Github link, and a methods writeup.
Submission limit: 6
Scoring
This is a foundational track — the goal is to validate that your pipeline can recover a clinically confirmed variant from real patient data. Prediction files are scored automatically on two primary metrics:

Rank points - How high did you rank the true variant(s). If your top-ranked variant(s) is correct, you'll get 100pts. If you get the correct variant(s) within your top 10 ranked ones, you'll get partial credit. If you correctly identify one of the two variants in a heterozygous pair, you'll get half credit.
F-max - the best balance of precision and recall your submission achieves at any confidence threshold, rewarding submissions that pinpoint the right answer without burying it in noise. Secondary or incidental findings won't hurt your automated score.
Your methods write-up is reviewed separately by a judging panel for scientific rigor and depth of understanding.
