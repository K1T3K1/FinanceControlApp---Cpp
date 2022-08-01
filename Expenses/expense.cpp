#include "expense.h"
using namespace Expenses;

Expense::Expense(float_precision cost, tm *date, EPaymentMethod pmethod)
{
    this->expenseCost = cost;
    this->expenseDate = date;
    this->paymentMethod = pmethod;
}

float_precision Expense::getExpenseCost()
{
    return this->expenseCost;
}

void Expense::setExpenseCost(float_precision costToSet)
{
    this->expenseCost = costToSet;
}

void Expense::getExpenseDate(tm *timeptr)
{
    timeptr = this->expenseDate;
}
void Expense::setExpenseDate(tm *timeptr)
{
    this->expenseDate = timeptr;
}
EPaymentMethod Expense::getPaymentMethod()
{
    return this->paymentMethod;
}