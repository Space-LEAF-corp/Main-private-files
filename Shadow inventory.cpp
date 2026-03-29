// ShadowInventory.cpp
#include "ShadowInventory.h"
#include "IItem.h"

bool UShadowInventory::AddItem_Implementation(TScriptInterface<IItem> Item, int32 Quantity)
{
    if (!Item)
        return false;

    FName Id = IItem::Execute_GetItemId(Item.GetObject());
    int64& Current = Storage.FindOrAdd(Id);
    Current += Quantity;
    return true;
}

bool UShadowInventory::RemoveItem_Implementation(FName ItemId, int32 Quantity)
{
    int64* Current = Storage.Find(ItemId);
    if (!Current || *Current < Quantity)
        return false;

    *Current -= Quantity;
    return true;
}

int32 UShadowInventory::GetQuantity_Implementation(FName ItemId) const
{
    const int64* Current = Storage.Find(ItemId);
    return Current ? static_cast<int32>(*Current) : 0;
}
