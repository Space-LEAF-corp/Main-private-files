// GameInventory.cpp
#include "GameInventory.h"
#include "IItem.h"

bool UGameInventory::AddItem_Implementation(TScriptInterface<IItem> Item, int32 Quantity)
{
    if (!Item)
        return false;

    FName Id = IItem::Execute_GetItemId(Item.GetObject());
    int32 MaxStack = IItem::Execute_GetMaxStack(Item.GetObject());

    int32& Current = Items.FindOrAdd(Id);
    if (Current + Quantity > MaxStack)
        return false;

    Current += Quantity;
    return true;
}

bool UGameInventory::RemoveItem_Implementation(FName ItemId, int32 Quantity)
{
    int32* Current = Items.Find(ItemId);
    if (!Current || *Current < Quantity)
        return false;

    *Current -= Quantity;
    return true;
}

int32 UGameInventory::GetQuantity_Implementation(FName ItemId) const
{
    const int32* Current = Items.Find(ItemId);
    return Current ? *Current : 0;
}
