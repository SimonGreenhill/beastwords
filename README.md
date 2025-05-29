# BEASTWORDS

Convert single partition beast2 linguistic XML files into word-partitioned analyses.


## Installation

Clone repository and run pip:

```shell
git clone https://github.com/SimonGreenhill/beastwords
cd beastwords && python -m pip install -e .
```

## Usage:


### Split on words

```shell
beastwords covarion.xml covarion.words.xml
```

### Divide into (n) partitions:

```shell
beastwords -p 5 covarion.xml covarion.5parts.xml
```

### Divide into specified groups:

Use this to split into groups defined by site position
This will divide the dataset into two partitions:
 p1 = 1-10
 p2 = 11-20

```shell
beastwords -p 1-10,11-20 covarion.xml covarion.groups.xml
```


## beastsitedistr can help you choose sizes:

Print a histogram of current partition sizes. In the below figure, there are 11 words with 11 sites (=cognate sets). 

```shell
beastsitedistr myfile.xml 

2	11	███████████
3	11	███████████
4	23	███████████████████████
5	17	█████████████████
6	19	███████████████████
7	17	█████████████████
8	18	██████████████████
9	17	█████████████████
10	15	███████████████
11	14	██████████████
12	10	██████████
13	9	█████████
14	6	██████
15	2	██
16	3	███
17	3	███
18	2	██
19	3	███
20	1	█
21	0
22	0
23	1	█
24	1	█
```
