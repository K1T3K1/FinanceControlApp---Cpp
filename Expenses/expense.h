#ifndef EXPENSE_H_
#define EXPENSE_H_
#include <stdlib.h>
#include <stddef.h>
#include "fprecision.h"

typedef struct tm tm;

namespace Expenses
{
    class Expense
    {
    public:
        Expense(float_precision cost, tm *date, EPaymentMethod pmethod);
        float_precision getExpenseCost(void);
        void setExpenseCost(float_precision costToSet);
        void getExpenseDate(tm *timeptr);
        void setExpenseDate(tm *timeptr);
        EPaymentMethod getPaymentMethod(void);
        void setPaymentMethod(EPaymentMethod methodToSet);

    protected:
        // BMP expenseTypeImage;
        float_precision expenseCost;
        tm *expenseDate;
        EPaymentMethod paymentMethod;
    };

    typedef enum
    {
        cash,
        card,
    } EPaymentMethod;

#endif // EXPENSE_H_
}