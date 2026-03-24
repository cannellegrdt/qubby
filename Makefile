##
## Project: qubby
## File name: Makefile
## Author: Cannelle Gourdet - lankley
## File description: Build system for the qubby quantum circuit simulator
##

NAME	:=	qubby

CXX	:=	g++
CXXFLAGS	:=	-std=c++20 -Wall -Wextra -Wpedantic
OPTFLAGS	:=	-O2
OMPFLAGS	:=	-fopenmp
DEBUGFLAGS	:=	-g3 -fsanitize=address,undefined

SRC_DIR	:=	src
INC_DIR	:=	include
OBJ_DIR	:=	obj
TEST_DIR	:=	tests
UNIT_DIR	:=	$(TEST_DIR)/unit
FUNC_DIR	:=	$(TEST_DIR)/functional

SRCS	:=	$(shell find $(SRC_DIR) -name '*.cpp' ! -name '*.sycl.cpp' 2>/dev/null)
OBJS	:=	$(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

SYCL_ALL_SRCS	:=	$(shell find $(SRC_DIR) -name '*.cpp' 2>/dev/null)

UNIT_BIN	:=	$(TEST_DIR)/unit_tests
FUNC_SCRIPT	:=	$(FUNC_DIR)/functional_tests.sh

UNIT_SRCS	:=	$(shell find $(UNIT_DIR) -name '*.cpp' 2>/dev/null)

LIBS	:=	-lcriterion
INC_FLAGS	:=	-I$(INC_DIR)

COVFLAGS	:=	--coverage -fprofile-arcs -ftest-coverage
COV_DIR	:=	coverage

SYCL_CXX	:=	icpx
SYCL_FLAGS	:=	-fsycl -std=c++20 -Wall

SYCL_OBJ_DIR	:=	obj/sycl
SYCL_TEST_BIN	:=	$(TEST_DIR)/sycl_tests
SYCL_SRCS	:=	$(filter-out $(SRC_DIR)/main.cpp, $(SYCL_ALL_SRCS))
SYCL_OBJS	:=	$(patsubst $(SRC_DIR)/%.cpp, $(SYCL_OBJ_DIR)/%.o, $(SYCL_SRCS))


all:	$(NAME)

$(NAME):	$(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(OMPFLAGS) $^ -o $@

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(OMPFLAGS) $(INC_FLAGS) -c $< -o $@

re:	fclean	all

clean:
	$(RM) -r $(OBJ_DIR)
	$(RM) -r $(COV_DIR)
	find . \( -name '*.gcno' -o -name '*.gcda' -o -name '*.gcov' \) | xargs $(RM) -f 2>/dev/null; true

fclean:	clean
	$(RM) $(NAME) $(UNIT_BIN) $(SYCL_TEST_BIN)
	$(RM) -r $(SYCL_OBJ_DIR)

unit_tests:	$(UNIT_BIN)

$(UNIT_BIN):	$(UNIT_SRCS) $(filter-out $(OBJ_DIR)/main.o, $(OBJS))
	@mkdir -p $(TEST_DIR)
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(INC_FLAGS) $^ $(LIBS) -o $@

functional_tests:	all
	@bash $(FUNC_SCRIPT)

run_tests: unit_tests
	@echo "── Unit tests ──────────────────────────────────────────"
	./$(UNIT_BIN) --color=always
	@echo "\n"
	@bash $(FUNC_SCRIPT)

coverage:	CXXFLAGS += $(COVFLAGS)
coverage:	OPTFLAGS  = -O0
coverage:	fclean $(UNIT_BIN)
	./$(UNIT_BIN)
	@mkdir -p $(COV_DIR)
	lcov --capture --directory . --output-file $(COV_DIR)/coverage.info \
	     --exclude '/usr/*' --exclude '$(TEST_DIR)/*'
	genhtml $(COV_DIR)/coverage.info --output-directory $(COV_DIR)/html
	@echo "\n Coverage report: $(COV_DIR)/html/index.html"

debug:	CXXFLAGS += $(DEBUGFLAGS)
debug:	OPTFLAGS  = -O0
debug:	re

sycl:
	@command -v $(SYCL_CXX) >/dev/null 2>&1 || \
	  { echo "Error: $(SYCL_CXX) not found — install Intel oneAPI Base Toolkit first."; exit 1; }
	$(SYCL_CXX) $(SYCL_FLAGS) $(INC_FLAGS) -DSYCL_ENABLED $(SYCL_ALL_SRCS) -o $(NAME)_sycl

$(SYCL_OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(SYCL_CXX) $(SYCL_FLAGS) $(INC_FLAGS) -DSYCL_ENABLED -c $< -o $@

$(SYCL_TEST_BIN):	$(UNIT_DIR)/test_sycl.cpp $(SYCL_OBJS)
	@mkdir -p $(TEST_DIR)
	$(SYCL_CXX) $(SYCL_FLAGS) $(INC_FLAGS) -DSYCL_ENABLED $^ $(LIBS) -o $@

sycl_tests:
	@command -v $(SYCL_CXX) >/dev/null 2>&1 || \
	  { echo "Error: $(SYCL_CXX) not found — install Intel oneAPI Base Toolkit first."; exit 1; }
	@$(MAKE) --no-print-directory $(SYCL_TEST_BIN)
	./$(SYCL_TEST_BIN) --color=always

format:
	find $(SRC_DIR) $(INC_DIR) $(TEST_DIR) -name '*.cpp' -o -name '*.hpp' \
	  | xargs clang-format -i --style=file 2>/dev/null || \
	  echo "clang-format not found — skipping"

lint:
	clang-tidy $(SRCS) -- $(CXXFLAGS) $(INC_FLAGS) 2>/dev/null || \
	  echo "clang-tidy not found — skipping"

help:
	@echo "Targets:"
	@echo "  all              Build $(NAME)"
	@echo "  re               Clean rebuild"
	@echo "  clean            Remove objects & coverage data"
	@echo "  fclean           clean + remove binaries"
	@echo "  unit_tests       Build unit test binary"
	@echo "  functional_tests Build functional test binary"
	@echo "  run_tests        Build & run both test suites"
	@echo "  coverage         Build with gcov, run tests, generate HTML report"
	@echo "  debug            Build with ASan/UBSan + no optimisations"
	@echo "  sycl             Build with Intel oneAPI (icpx -fsycl)"
	@echo "  sycl_tests       Build & run SYCL unit tests (requires icpx)"
	@echo "  format           Run clang-format on all sources"
	@echo "  lint             Run clang-tidy on all sources"


.PHONY: all re clean fclean	\
		unit_tests functional_tests run_tests coverage	\
		debug sycl sycl_tests format lint help
