import torch
from scipy.stats import mode

from .utils import ancestry_accuracy, ProgressSaver, AverageMeter, ReshapedCrossEntropyLoss,\
    adjust_learning_rate, to_device, correct_max_indices, compute_ibd

import time

# def train(model, train_loader, valid_loader, args):

#     device  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
#     print("device:", device)

#     model.to(device)

#     criterion = ReshapedCrossEntropyLoss()


#     # Basic
#     optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

#     # Different LR for topK weights
#     # optimizer = torch.optim.Adam([
#     #     {'params':model.smoother.parameters()},
#     #     {'params':model.add_poolings.parameters(), 'lr':args.lr/10}
#     #     ], lr=args.lr)
#     # optimizer = torch.optim.SGD(model.parameters(), lr=args.lr, momentum=0.9)

#     init_time = time.time()

#     progress_saver = ProgressSaver(args.exp)
#     train_loss_meter = AverageMeter()
#     best_val_acc = -1
#     best_epoch = -1

#     lr = args.lr


#     init_epoch = 0
#     if args.resume:
#         progress_saver.load_progress()
#         init_epoch, best_val_loss, start_time = progress_saver.get_resume_stats()

#         init_time = time.time() - start_time
#         model.load_state_dict(torch.load(args.exp + "/models/last_model.pth"))
#         optimizer.load_state_dict(
#             torch.load(args.exp + "/models/last_optim.pth"))
#         for state in optimizer.state.values():
#             for k, v in state.items():
#                 if isinstance(v, torch.Tensor):
#                     state[k] = v.to(device)
#         print("loaded state dict from epoch %d" % init_epoch)

#         init_epoch += 1

#     model.to(device)
#     for n in range(init_epoch, args.num_epochs):
#         model.train()
#         train_loss_meter.reset()

#         if args.lr_decay > 0:
#             lr = adjust_learning_rate(args.lr, args.lr_decay, optimizer, n)

#         for i, batch in enumerate(train_loader):

#             batch = to_device(batch, device)

#             output = model(batch["mixed_vcf"], batch["ref_panel"])
#             loss = criterion(output["predictions"], batch["mixed_labels"].to(device))
#             loss.backward()
#             if((i+1) % args.update_every) == 0:
#                 optimizer.step()
#                 optimizer.zero_grad()

#             train_loss_meter.update(loss.item())

#         val_acc, val_loss = validate(model, valid_loader, criterion, args)
#         train_loss = train_loss_meter.get_average()

#         total_time = time.time() - init_time

#         if best_val_acc < val_acc:
#             best_val_acc = val_acc
#             best_val_loss = val_loss
#             best_epoch = n
#             torch.save(model.state_dict(), args.exp + "/models/best_model.pth")

#         torch.save(model.state_dict(), args.exp + "/models/last_model.pth")
#         torch.save(optimizer.state_dict(), args.exp + "/models/last_optim.pth")

#         epoch_data = {
#             "epoch": n,
#             "train_loss": train_loss,
#             "val_loss": val_loss,
#             "val_acc": val_acc.cpu(),
#             "best_epoch": best_epoch,
#             "best_val_loss": best_val_loss,
#             "time": total_time,
#             "lr": lr
#         }

#         progress_saver.update_epoch_progess(epoch_data)

#         print("epoch #", n, ":\tVal acc:", val_acc.item(), "\ttime:", time.time()- init_time)

def train(model, train_loader, valid_loader, args):

    device  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", device)

    model.to(device)

    criterion = ReshapedCrossEntropyLoss()

    # Basic optimizer
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    init_time = time.time()

    progress_saver = ProgressSaver(args.exp)
    train_loss_meter = AverageMeter()
    best_val_acc = -1
    best_epoch = -1

    lr = args.lr

    init_epoch = 0
    if args.resume:
        progress_saver.load_progress()
        init_epoch, best_val_loss, start_time = progress_saver.get_resume_stats()

        init_time = time.time() - start_time
        #model.load_state_dict(torch.load(args.exp + "/models/last_model.pth"))
        state_dict = torch.load(args.exp + "/models/best_model.pth") 
        model_state_dict = model.state_dict()
        filtered_state_dict = {k: v for k, v in state_dict.items() if k in model_state_dict}

        model_state_dict.update(filtered_state_dict)  # Update the model's state_dict with the loaded weights
        model.load_state_dict(model_state_dict)  # Load the updated state_dict into the model

        # model.load_state_dict(torch.load(args.exp + "/models/best_model.pth"))
        optim = torch.load(args.exp + "/models/last_optim.pth")

        try:
            # Load the optimizer's state dictionary if possible
            optimizer.load_state_dict(optim)
            print("Optimizer state loaded successfully.")
        except ValueError as e:
            print(f"Error loading optimizer state: {e}")
            # Handle this error by reinitializing the optimizer from scratch
            optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
            print("Reinitialized optimizer with new parameters.")

        # Move optimizer states to the correct device (if necessary)
        for state in optimizer.state.values():
            for k, v in state.items():
                if isinstance(v, torch.Tensor):
                    state[k] = v.to(device)
        print("loaded state dict from epoch %d" % init_epoch)

        init_epoch += 1

    model.to(device)

    for n in range(init_epoch, args.num_epochs):
        model.train()
        train_loss_meter.reset()

        # Track correct predictions for accuracy
        correct_predictions = 0
        total_samples = 0

        if args.lr_decay > 0:
            lr = adjust_learning_rate(args.lr, args.lr_decay, optimizer, n)

        for i, batch in enumerate(train_loader):
            batch = to_device(batch, device)

            output = model(batch["mixed_vcf"], batch["ref_panel"])
            snp_predictions = output["predictions"]  # Shape: [batch_size, num_snps, num_classes]

            # Aggregating predictions across SNPs (mode along SNP axis)
            snp_predictions = snp_predictions.argmax(dim=-1)  # Shape: [batch_size, num_snps]
            mode_predictions, _ = mode(snp_predictions.cpu().numpy(), axis=1)  # Mode along SNP axis

            # Aggregating mixed labels across SNPs (mode along SNP axis)
            mixed_labels = batch["mixed_labels"]  # Shape: [batch_size, num_snps]
            mode_labels, _ = mode(mixed_labels.cpu().numpy(), axis=1)  # Mode across SNPs

            # Convert mode predictions and labels back to tensor
            mode_predictions = torch.tensor(mode_predictions, dtype=torch.long, device=device)
            mode_labels = torch.tensor(mode_labels, dtype=torch.long, device=device)

            # Calculate correct predictions
            correct_predictions += (mode_predictions == mode_labels).sum().item()
            total_samples += mode_labels.size(0)

            # Compute loss
            loss = criterion(output["predictions"], batch["mixed_labels"].to(device))
            loss.backward()
            if((i+1) % args.update_every) == 0:
                optimizer.step()
                optimizer.zero_grad()

            train_loss_meter.update(loss.item())

        # Calculate training accuracy
        train_accuracy = correct_predictions / total_samples
        print("train_accuracy",train_accuracy)

        # Validate the model after each epoch
        val_acc, val_loss = validate(model, valid_loader, criterion, args)
        train_loss = train_loss_meter.get_average()

        total_time = time.time() - init_time

        # Save the model if the validation accuracy improves
        if best_val_acc < val_acc:
            best_val_acc = val_acc
            best_val_loss = val_loss
            best_epoch = n
            torch.save(model.state_dict(), args.exp + "/models/best_model.pth")

        torch.save(model.state_dict(), args.exp + "/models/last_model.pth")
        torch.save(optimizer.state_dict(), args.exp + "/models/last_optim.pth")

        epoch_data = {
            "epoch": n,
            "train_loss": train_loss,
            "val_loss": val_loss,
            "val_acc": val_acc,
            "best_epoch": best_epoch,
            "best_val_loss": best_val_loss,
            "time": total_time,
            "lr": lr
        }

        print(epoch_data)

        progress_saver.update_epoch_progess(epoch_data)

        print(f"epoch #{n}: Train accuracy: {train_accuracy:.4f}, Val acc: {val_acc:.4f}, time: {time.time() - init_time:.2f}s")


# def validate(model, val_loader, criterion, args):

#     with torch.no_grad():

#         val_loss = AverageMeter()

#         device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

#         model.eval().to(device)

#         acc = torch.tensor(0).float()
#         for i, batch in enumerate(val_loader):
#             batch = to_device(batch, device)

#             output = model(batch["mixed_vcf"], batch["ref_panel"])

#             acc = acc + ancestry_accuracy(output["predictions"], batch["mixed_labels"])
#             loss = criterion(output["predictions"], batch["mixed_labels"])
#             val_loss.update(loss.item())

#         acc = acc / len(val_loader.dataset)

#         return acc, val_loss.get_average()

def validate(model, val_loader, criterion, args):

    with torch.no_grad():
        val_loss = AverageMeter()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model.eval().to(device)

        correct_predictions = 0  # Initialize correct predictions count
        total_samples = 0  # Initialize total sample coun

        acc = torch.tensor(0).float()

        for i, batch in enumerate(val_loader):
            batch = to_device(batch, device)

            # Get the model predictions
            output = model(batch["mixed_vcf"], batch["ref_panel"])
            snp_predictions = output["predictions"]  # Shape: [batch_size, num_snps, num_classes]

            # Aggregating predictions across SNPs (mode along SNP axis)
            snp_predictions = snp_predictions.argmax(dim=-1)  # Shape: [batch_size, num_snps]
            mode_predictions, _ = mode(snp_predictions.cpu().numpy(), axis=1)  # Mode along SNP axis

            # Aggregating mixed labels across SNPs (mode along SNP axis)
            mixed_labels = batch["mixed_labels"]  # Shape: [batch_size, num_snps]
            mode_labels, _ = mode(mixed_labels.cpu().numpy(), axis=1)  # Mode across SNPs

            # Convert mode predictions and labels back to tensor
            mode_predictions = torch.tensor(mode_predictions, dtype=torch.long, device=device)
            mode_labels = torch.tensor(mode_labels, dtype=torch.long, device=device)

            # Calculate correct predictions

            correct_predictions += (mode_predictions == mode_labels).sum().item()
            total_samples += mode_labels.size(0)

                    # Calculate accuracy over all batches
            accuracy = correct_predictions / total_samples

            # Compute loss
            loss = criterion(output["predictions"], batch["mixed_labels"])
            val_loss.update(loss.item())


        # # Calculate accuracy over all batches
        # acc = acc / len(val_loader.dataset)

        return accuracy, val_loss.get_average()

def inference(model, test_loader, args):

    with torch.no_grad():

        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        model.eval().to(device)


        all_predictions = []
        all_predictions_window = []
        all_ibd = []

        for i, batch in enumerate(test_loader):
            batch = to_device(batch, device)

            output = model(batch["mixed_vcf"], batch["ref_panel"])

            output['max_indices'] = correct_max_indices(output['max_indices'], batch['reference_idx'])

            ibd = compute_ibd(output)

            predicted_labels = torch.argmax(output['predictions'], dim=2)
            predicted_labels_window = torch.argmax(output['out_smoother'], dim=1)

            all_predictions.append(predicted_labels)
            all_predictions_window.append(predicted_labels_window)
            all_ibd.append(ibd)

        all_predictions = torch.cat(all_predictions, dim=0)
        all_predictions_window = torch.cat(all_predictions_window, dim=0)
        all_ibd = torch.cat(all_ibd, dim=0)

        return all_predictions, all_predictions_window, all_ibd

