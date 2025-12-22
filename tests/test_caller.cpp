#include <gtest/gtest.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sys/wait.h>
#include <unistd.h>

// Placeholder test to verify Google Test infrastructure works
TEST(InfrastructureTest, PlaceholderTest) {
    EXPECT_TRUE(true);
}

// Test will be expanded in subsequent issues
TEST(InfrastructureTest, BasicAssertions) {
    EXPECT_EQ(1 + 1, 2);
    EXPECT_NE(1, 2);
}

// Test file open error handling
TEST(FileHandlingTest, NonExistentFileReturnsError) {
    std::string nonexistent = "/tmp/nonexistent_" + std::to_string(getpid()) + "_" + std::to_string(time(nullptr)) + ".fa";
    std::string cmd = "./bin/callerpp -i " + nonexistent + " 2>/dev/null";
    int result = std::system(cmd.c_str());
    ASSERT_NE(result, -1) << "std::system() failed to execute";
    ASSERT_TRUE(WIFEXITED(result)) << "Process did not exit normally";
    EXPECT_NE(WEXITSTATUS(result), 0);
}

TEST(FileHandlingTest, ValidFileWorks) {
    // Use unique filename with process ID to avoid race conditions
    std::string temp_path = "/tmp/test_callerpp_" + std::to_string(getpid()) + ".fa";

    // RAII cleanup guard ensures file is removed even if assertions fail
    struct FileGuard {
        std::string path;
        ~FileGuard() { std::remove(path.c_str()); }
    } guard{temp_path};

    // Create a temporary test file
    std::ofstream testfile(temp_path);
    ASSERT_TRUE(testfile.is_open()) << "Failed to create temporary test file";
    testfile << ">test\nACGT\nACGT\n";
    testfile.close();

    std::string cmd = "./bin/callerpp -i " + temp_path + " >/dev/null 2>&1";
    int result = std::system(cmd.c_str());
    ASSERT_NE(result, -1) << "std::system() failed to execute";
    ASSERT_TRUE(WIFEXITED(result)) << "Process did not exit normally";
    EXPECT_EQ(WEXITSTATUS(result), 0);
}

// Test input validation for numeric arguments
TEST(InputValidationTest, InvalidAlgorithmReturnsError) {
    // Algorithm must be 0, 1, or 2
    int result = std::system("echo '>test\nACGT' | ./bin/callerpp -a 5 2>/dev/null");
    EXPECT_NE(result, 0);
}

TEST(InputValidationTest, NonNumericAlgorithmReturnsError) {
    int result = std::system("echo '>test\nACGT' | ./bin/callerpp -a abc 2>/dev/null");
    EXPECT_NE(result, 0);
}

TEST(InputValidationTest, InvalidResortReturnsError) {
    // Resort must be 0, 1, or 2
    int result = std::system("echo '>test\nACGT' | ./bin/callerpp -r 10 2>/dev/null");
    EXPECT_NE(result, 0);
}

TEST(InputValidationTest, ValidAlgorithmWorks) {
    int result = std::system("echo '>test\nACGT\nACGT' | ./bin/callerpp -a 0 >/dev/null 2>&1");
    EXPECT_EQ(result, 0);
    result = std::system("echo '>test\nACGT\nACGT' | ./bin/callerpp -a 1 >/dev/null 2>&1");
    EXPECT_EQ(result, 0);
    result = std::system("echo '>test\nACGT\nACGT' | ./bin/callerpp -a 2 >/dev/null 2>&1");
    EXPECT_EQ(result, 0);
}
