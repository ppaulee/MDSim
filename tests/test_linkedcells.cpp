#include "gtest/gtest.h"

// Tests correct construction of cell grid
TEST(LinkedCells, gridConstruction) {
    auto cells = new LinkedCells({2,2,2}, 1, 1);
    EXPECT_EQ(1, 1);
    for (int x = 0; x < 6; ++x) {
        for (int y = 0; y < 6; ++y) {
            for (int z = 0; z < 6; ++z) {
                bool halo = (x == 0 || y == 0 || z == 0 || x == 5 || y == 5 || z == 5);
                bool boundary = (x == 1 || y == 1 || z == 1 || x == 4 || y == 4 || z == 4);

                EXPECT_EQ(halo, cells->isHaloCell({x,y,z}));
                EXPECT_EQ(boundary, cells->isBoundaryCell({x,y,z}));
            }
        }
    }
}