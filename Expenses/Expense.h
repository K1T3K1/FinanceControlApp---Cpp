#ifndef EXPENSE_H_
#define EXPENSE_H_
#include <stdlib.h>
#include <stddef.h>
#include <string>
#include "fprecision.h"

typedef struct tm tm;

namespace Expenses
{
    class Expense
    {
    public:
        Expense(std::string title, float_precision cost, tm *date);
        void setExpenseTitle(std::string titleToSet);
        std::string getExpenseTitle(void);
        float_precision getExpenseCost(void);
        void setExpenseCost(float_precision costToSet);
        void getExpenseDate(tm *timeptr);
        void setExpenseDate(tm *timeptr);

    protected:
        // BMP expenseTypeImage;
        std::string expenseTitle;
        float_precision expenseCost;
        tm *expenseDate;
    };
}

#endif // EXPENSE_H_