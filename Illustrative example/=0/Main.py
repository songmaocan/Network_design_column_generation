"""
Algorithm: Solve the reliability-oriented network design problem by Column generation
Copyright: Maocan Song, 1097133316@qq.com
"""
import time
from Method import Solve

def main():
    start_time = time.time()
    mod=Solve()
    mod.g_solving_RNDP_by_CG()
    end_time = time.time()
    spend_time = end_time - start_time
    mod.output_results(spend_time)
    print("CPU running time {} min".format(round(spend_time / 60), 3))

if __name__ == '__main__':
    main()