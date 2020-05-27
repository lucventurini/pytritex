import numpy as np
from .break_into_blocks import break_into_blocks
import array as carray


def block_optimise(edges, path):
    block_chain = break_into_blocks(edges, path)
    if len(block_chain.blocks) == path.shape[0]:
        return path
    elif len(block_chain.blocks) == 1:  # Only one block, no improvement possible
        return path

    block_distance = carray.array("i")