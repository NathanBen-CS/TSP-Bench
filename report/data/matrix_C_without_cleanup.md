| Method | uniform 1k | uniform 10k | uniform 100k | uniform 1M | clustered 1k | clustered 10k | clustered 100k |
|---|---|---|---|---|---|---|---|
| nearest neighbour (nn) | 25.5 | 24.8 | 23.6 | - | 28.9 | 30.1 | 32.0 |
| greedy (DIMACS greedy) | 16.6 | 16.4 | 14.7 | 14.3 | 18.0 | 23.1 | 21.7 |
| greedy start of our solver | 17.2 | 16.4 | 15.1 | 14.3 | 17.3 | 21.1 | 21.1 |
| Hilbert order (spacefill) | 38.3 | 38.6 | 38.5 | 38.2 | 50.1 | 55.6 | 55.9 |
| Moore order | 34.8 | 37.1 | 38.0 | 38.1 | 42.4 | 52.3 | 54.9 |
| Sierpinski order | 31.5 | 34.5 | 34.9 | 35.1 | 38.5 | 50.0 | 52.1 |
| Knopp order | - | 35.6 | - | - | - | 52.7 | - |
| Peano order | - | 47.5 | - | - | - | 66.4 | - |
| Morton order | - | 103.1 | - | - | - | 113.2 | - |
| Strip order | - | 4,584 | - | - | - | 7,321 | - |
| hull rings in order | 334.3 | 562.4 | 885.9 | 1,353 | 884.2 | 1,718 | 2,977 |
| boxes + merge, hilbert | 16.3 (M4) | 15.9 (M4) | 16.9 (M4) | 16.8 (M4) | 28.4 (M4) | 28.8 (M5) | 29.6 (M5) |
| boxes + merge, moore | 16.3 (M4) | 15.9 (M4) | 16.9 (M4) | 16.8 (M4) | 24.3 (M4) | 29.3 (M4) | 29.7 (M5) |
| boxes + merge, sierpinski | 20.0 (M4) | 20.2 (M4) | 19.9 (M4) | 20.1 (M4) | 28.9 (M5) | 32.2 (M4) | 31.8 (M4) |
| boxes + merge, knopp | 20.0 (M4) | 20.2 (M4) | 19.9 (M4) | 20.1 (M4) | 24.5 (M4) | 32.6 (M5) | 32.0 (M5) |
| boxes + merge, peano | 15.6 (M4) | 16.9 (M4) | 17.2 (M4) | 17.4 (M4) | 24.2 (M4) | 30.2 (M5) | 31.0 (M5) |
| boxes + merge, morton | 16.3 (M4) | 15.9 (M4) | 16.9 (M4) | 16.8 (M4) | 26.8 (M4) | 32.1 (M5) | 31.6 (M5) |
| boxes + merge, strip | 284.3 (M4) | 1,014 (M4) | 3,409 (M4) | 28,502 (M1) | 146.2 (M4) | 1,046 (M4) | 5,487 (M4) |
| hull rings + merge | 269.4 (M4) | 505.4 (M4) | 840.1 (M4) | 1,318 (M4) | 451.3 (M4) | 1,255 (M4) | 2,547 (M4) |
