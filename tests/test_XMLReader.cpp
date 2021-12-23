//
// Created by ethan on 12/23/2021.
//

#include "gtest/gtest.h"
#include "../src/XMLReader/driver.h"
#include "../src/XMLReader/library.h"

TEST(driver, generateFromXML){

    library::molsim m = generateFromXML("../tests/testinput.xml");
    EXPECT_EQ(m.input(), "../input.txt");
    EXPECT_EQ(m.delta_t(), 1);
    EXPECT_EQ(m.endtime(), 2);
    EXPECT_EQ(m.outputStep(), 3);
    EXPECT_EQ(m.epsilon(), 4);
    EXPECT_EQ(m.sigma(), 5);

}