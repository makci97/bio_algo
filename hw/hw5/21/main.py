import sys
import numpy as np
from functools import lru_cache
from collections import Counter, defaultdict


def min_coins_num(money, coins):
    # init dp
    dp = defaultdict(lambda: np.inf)
    dp[0] = 0
    for i in range(1, money + 1):
        dp[i] = 1 + min(dp[i - coin] for coin in coins)

    return dp[money]


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        money = int(file.readline().strip())
        coins = tuple(map(int, file.readline().strip().split(',')))

    # write
    with open(output_filename, 'w') as file:
        file.write(str(min_coins_num(money, coins)))
