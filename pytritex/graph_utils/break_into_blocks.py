import numpy as np


class Block:
    def __init__(self, orientation, markers, size, first, last, previous_block_index,
                 next_block_index):
        self.orientation = orientation, self.markers = orientation, markers
        self.size, self.first, self.last = size, first, last
        self.previous_block_index, self.next_block_index = previous_block_index, next_block_index


class BlockChain:
    def __init__(self):
        self.header = 0
        self.blocks = []


def break_into_blocks(edge_list, path: np.array):

    """Reimplementation of MSTOpt::break_into_blocks from MstMap"""
    # TODO check the implementation of what the BlockChain is
    block_chain = BlockChain()
    last_block_break = 0
    for position in np.arange(1, path.shape[0] - 1, dtype=np.int):
        before, node, after = path[[position - 1, position, position + 1]]
        cost_of_edge = edge_list[[before, node]]
        nodes_before = path[np.arange(0, position)]
        nodes_after = path[np.arange(position + 1, path.shape[0])]
        valid = True
        # TODO: maybe it's possible to vectorise this
        for before_node in nodes_before:
            if edge_list[[before_node, node]] < cost_of_edge:
                valid = False
                break
        for after_node in nodes_after:
            if edge_list[[node, after_node]] < cost_of_edge:
                valid = False
                break
        if valid is False:
            # Construct a new block
            # Block crt_block = {true, block_markers, block_markers.size(),
            #              block_markers[0], block_markers[block_markers.size() - 1],
            #              -1, -1};
            # Save that the last position we have broken up is here
            block_markers = path[np.arange(last_block_break, position - 1, dtype=np.int)]
            current_block = Block(True, block_markers, block_markers.shape[0],
                                  block_markers[0], block_markers[-1], -1, -1)
            block_chain.blocks.append(current_block)
            last_block_break = position
    # Assemble the last block
    last_block_markers = path[np.arange(last_block_break, path.shape[0], dtype=np.int)]
    last_block = Block(True, last_block_markers, last_block_markers.shape[0],
                       last_block_markers[0], last_block_markers[-1], -1, -1)
    block_chain.blocks.append(last_block)
    # Do some correctness checks
    total_bins = 0
    for position, block in enumerate(block_chain.blocks):
        total_bins += block.markers.shape[0]
        block_chain.blocks[position].previous_block_index = position - 1
        block_chain.blocks[position].next_block_index = position + 1
    assert total_bins == path.shape[0]
    block_chain.blocks[-1].next_block_index = -1
    block_chain.header = 0
    return block_chain
