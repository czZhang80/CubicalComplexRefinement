#include "test.h"
#include <exception>
#include <iostream>

int main(void)
{
    try
    {
        testMoveLeft();
        testMoveRight();
        testMoveUp();
        testMoveDown();
    }
    catch(const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}
